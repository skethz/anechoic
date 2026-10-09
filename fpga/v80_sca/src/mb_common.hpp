// Multi-bit coupling helpers (MULTIBIT_SPEC.md) shared by the vector generator (v6/gen_vectors_mb.cpp) and the board host
// (host_sca_multi_mb.cpp). The golden model stays sca_ref.hpp run_trial(..., dense) with its arithmetic unchanged: these
// helpers only build the dense int8 matrix, the bit planes streamed to the engine, the scaled tables and an independent
// scorer from the raw edges.
#pragma once
#include "sca_ref.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mb {

struct Edge { int i, j, w; };   // 0-based, i < j, raw Max-Cut weight w (Ising J_ij = -w)

struct Matrix {
    std::vector<int8_t> J = std::vector<int8_t>(size_t(sca::NP) * sca::NP, 0);   // J[x * NP + y]; symmetric, J_xx = 0, padding 0
    int64_t sumw = 0;          // sum_{i<j} w_ij with w = -J (cut = (sumw + sum_i s_i h_i / 2) / 2)
    int maxabs = 0;            // max |J|
    std::string source;        // description
    std::vector<Edge> edges;   // raw edges when loaded from an edge list (independent scorer); empty otherwise
    bool from_packed = false;  // K2000: the packed +-1 graph is kept for the original independent scorer
    sca::Graph packed;
};

inline int8_t& at(Matrix& m, int x, int y) { return m.J[size_t(x) * sca::NP + y]; }
inline int8_t at(const Matrix& m, int x, int y) { return m.J[size_t(x) * sca::NP + y]; }

inline void finish(Matrix& m) {
    int64_t s = 0; int mx = 0;
    for (int x = 0; x < sca::NP; ++x)
        for (int y = 0; y < sca::NP; ++y) {
            const int v = at(m, x, y);
            if (v != at(m, y, x)) throw std::runtime_error("matrix not symmetric");
            if ((x == y || x >= sca::N || y >= sca::N) && v != 0) throw std::runtime_error("nonzero diagonal or padding");
            mx = std::max(mx, std::abs(v));
            if (x < y) s -= v;
        }
    m.sumw = s; m.maxabs = mx;
}

// K2000 (+-1, packed "SNOWGPU1" file): identical dense matrix to sca::dense_matrix(load_graph(path)).
inline Matrix from_k2000(const std::string& path) {
    Matrix m; m.packed = sca::load_graph(path); m.from_packed = true; m.source = "k2000:" + path;
    const std::vector<int8_t> d = sca::dense_matrix(m.packed);
    m.J = d; finish(m);
    if (m.sumw != m.packed.sumw) throw std::runtime_error("K2000 sumw mismatch");
    return m;
}

// G-set (Stanford format): first line "n e", then e lines "i j w" (1-based). n must equal the engine's N (2000).
inline Matrix from_gset(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    int n = 0; long e = 0;
    if (!(f >> n >> e) || n != sca::N) throw std::runtime_error("G-set instance must have N = 2000 for this engine: " + path);
    Matrix m; m.source = "gset:" + path;
    for (long k = 0; k < e; ++k) {
        int i, j, w;
        if (!(f >> i >> j >> w)) throw std::runtime_error("short G-set file");
        --i; --j;
        if (i < 0 || j < 0 || i >= n || j >= n || i == j || w == 0 || std::abs(w) > 127) throw std::runtime_error("bad G-set edge");
        if (i > j) std::swap(i, j);
        if (at(m, i, j) != 0) throw std::runtime_error("duplicate G-set edge");
        at(m, i, j) = int8_t(-w); at(m, j, i) = int8_t(-w);
        m.edges.push_back({i, j, w});
    }
    std::string rest; if (f >> rest) throw std::runtime_error("trailing data in G-set file");
    finish(m);
    return m;
}

// Random symmetric matrix with N = 2000: each pair i<j nonzero with probability `density`, |J| uniform in 1..maxabs, random
// sign. splitmix64 stream, so the matrix is a pure function of (seed, density, maxabs).
inline Matrix random_k(uint64_t seed, double density, int maxabs) {
    Matrix m; m.source = "rand:" + std::to_string(seed) + ":" + std::to_string(density) + ":" + std::to_string(maxabs);
    uint64_t st = seed;
    const uint64_t thr = uint64_t(density * 18446744073709551616.0 - 1.0);
    for (int i = 0; i < sca::N; ++i)
        for (int j = i + 1; j < sca::N; ++j) {
            const uint64_t a = sca::Rng::splitmix(st), b = sca::Rng::splitmix(st);
            if (a > thr) continue;
            const int mag = 1 + int(b % uint64_t(maxabs)); const int v = ((b >> 32) & 1) ? mag : -mag;
            at(m, i, j) = int8_t(v); at(m, j, i) = int8_t(v);
            m.edges.push_back({i, j, -v});
        }
    finish(m);
    return m;
}

// Graph carrying sumw for run_trial (its packed J is unused when the dense matrix is given).
inline sca::Graph graph_of(const Matrix& m) {
    sca::Graph g = m.from_packed ? m.packed : sca::Graph{};
    g.sumw = m.sumw;
    return g;
}

// Coupling planes for the K-bit engine, plane-major: word (p * NP + x) * 32 + w holds columns 64w..64w+63 of row x of plane
// p (bit b = column 64w + b). Plane 0 = sign (1 iff J < 0), plane p >= 1 = bit p-1 of |J|.
inline std::vector<uint64_t> planes(const Matrix& m, int K) {
    if (K < 2 || K > 8 || m.maxabs > (1 << (K - 1)) - 1) throw std::runtime_error("matrix does not fit K planes");
    std::vector<uint64_t> P(size_t(K) * sca::NP * sca::WORDS, 0);
    for (int x = 0; x < sca::NP; ++x)
        for (int y = 0; y < sca::NP; ++y) {
            const int v = at(m, x, y); if (v == 0) continue;
            const int a = std::abs(v);
            for (int p = 0; p < K; ++p) {
                const bool bit = p == 0 ? v < 0 : ((a >> (p - 1)) & 1) != 0;
                if (bit) P[(size_t(p) * sca::NP + x) * sca::WORDS + (y >> 6)] |= uint64_t(1) << (y & 63);
            }
        }
    return P;
}

// MULTIBIT_SPEC correction scale: gamma_t = lambda_t * (mean_i sum_j J_ij^2) / (N * 2 T_{t-1}); this returns
// (mean_i sum_j J_ij^2) / N. (+-1 complete graph: (N - 1) / N; the K2000 regression uses scale 1, i.e. v6.4's tables.)
inline double j2_scale(const Matrix& m) {
    double s = 0;
    for (int x = 0; x < sca::N; ++x)
        for (int y = 0; y < sca::N; ++y) { const double v = at(m, x, y); s += v * v; }
    return s / sca::N / sca::N;
}

// Schedule -> fixed-point tables (sca::make_tables, unchanged), with the correction scale folded into lambda (Onsager) or
// kappa (TEC-T, kconst = kappa * T_{t-1} * ramp in Q16.16, as host_sca_multi_v6.cpp). scale == 1 gives v6.4's tables.
struct Sched {
    double t0 = 30, t1 = 5, q = 4, lam = 0, jv = 0, kappa = 0, scale = 1;
    int S = 560; bool ramp = false;
};
inline std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S);
    for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1));
    return T;
}
inline sca::Tables tables(const Sched& c) {
    const std::vector<double> T = geom(c.t0, c.t1, c.S);
    std::vector<double> L(c.S, c.lam);
    if (c.ramp) for (int t = 0; t < c.S; ++t) if (t >= 0.7 * c.S) L[t] = c.lam * (c.S - 1 - t) / std::max(1.0, c.S - 1 - 0.7 * c.S);
    if (c.scale != 1.0) for (double& l : L) l *= c.scale;
    sca::Tables tb = sca::make_tables(T, std::vector<double>(c.S, c.q), L, c.jv);
    if (c.kappa > 0)
        for (int t = 1; t < c.S; ++t) {
            const double rf = (c.ramp && t >= 0.7 * c.S) ? (c.S - 1 - t) / std::max(1.0, c.S - 1 - 0.7 * c.S) : 1.0;
            tb.kconst[t] = int32_t(std::llround(c.kappa * c.scale * T[t - 1] * rf * 65536.0));
            tb.kcorr[t] = 0;
        }
    return tb;
}

// Engine range checks (sca_core_mb.v): kcorr fits int32; 4T fits the 27-bit lane constant; |q| + |corr| < 2^29 so that
// z = s*h*2^16 + C stays exact in the 32-bit (K <= 4) lane representation.
inline void check_tables(const sca::Tables& tb, int K) {
    for (size_t t = 0; t < tb.fourT.size(); ++t) {
        if (tb.kcorr[t] != int64_t(int32_t(tb.kcorr[t]))) throw std::runtime_error("kcorr does not fit int32");
        if (tb.fourT[t] <= 0 || tb.fourT[t] >= (1 << 27)) throw std::runtime_error("4T outside the 27-bit lane range");
        const int64_t corr = ((int64_t(2048) * std::llabs(tb.kcorr[t])) >> 8) + std::llabs(int64_t(tb.kconst[t]));
        if (K <= 4 && std::llabs(int64_t(tb.q[t])) + corr >= (int64_t(1) << 29)) throw std::runtime_error("|q| + |corr| >= 2^29");
    }
}

// Independent scorer: K2000 from the packed +-1 graph (as host_sca_multi_v6.cpp), otherwise from the raw edge list.
inline int64_t cut_of(const Matrix& m, const uint64_t* bits) {
    auto sp = [&](int i) { return ((bits[i >> 6] >> (i & 63)) & 1) != 0; };
    int64_t c = 0;
    if (m.from_packed) {
        for (int x = 0; x < sca::N; ++x)
            for (int y = x + 1; y < sca::N; ++y)
                if (sp(x) != sp(y)) c += -sca::jval(m.packed, x, y);
        return c;
    }
    if (m.edges.empty()) throw std::runtime_error("no raw edges for independent scoring");
    for (const Edge& e : m.edges) if (sp(e.i) != sp(e.j)) c += e.w;
    return c;
}

// FNV-1a 64 of the dense matrix (identity check in logs; file SHA-256 values are recorded by the scripts)
inline uint64_t fingerprint(const Matrix& m) {
    uint64_t h = 1469598103934665603ull;
    for (int8_t v : m.J) { h ^= uint8_t(v); h *= 1099511628211ull; }
    return h;
}

inline Matrix load(const std::string& spec) {
    if (spec.rfind("k2000:", 0) == 0) return from_k2000(spec.substr(6));
    if (spec.rfind("gset:", 0) == 0) return from_gset(spec.substr(5));
    if (spec.rfind("rand:", 0) == 0) {   // rand:<seed>:<density>:<maxabs>
        std::stringstream ss(spec.substr(5)); std::string a, b, c;
        std::getline(ss, a, ':'); std::getline(ss, b, ':'); std::getline(ss, c, ':');
        return random_k(std::stoull(a), std::stod(b), std::stoi(c));
    }
    throw std::runtime_error("matrix spec must be k2000:<file>, gset:<file> or rand:<seed>:<density>:<maxabs>");
}

}  // namespace mb
