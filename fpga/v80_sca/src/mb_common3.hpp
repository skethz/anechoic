// r3 helpers (copy of mb_common.hpp, which stays frozen for PROTOCOL_HW_MB): runtime active-spin count n, per-spin biases
// (MULTIBIT_SPEC.md "Bias (external field)"), a text Ising format and the bias stream segment. Golden model for r3:
// sca_ref_bias.hpp run_trial_bias (unchanged shared file). Original header text follows.
// Multi-bit coupling helpers (MULTIBIT_SPEC.md) shared by the vector generator (v6/gen_vectors_mb.cpp) and the board host
// (host_sca_multi_mb.cpp). The golden model stays sca_ref.hpp run_trial(..., dense) with its arithmetic unchanged: these
// helpers only build the dense int8 matrix, the bit planes streamed to the engine, the scaled tables and an independent
// scorer from the raw edges.
#pragma once
#include "sca_ref_bias.hpp"
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

namespace mb3 {

struct Edge { int i, j, w; };   // 0-based, i < j, raw Max-Cut weight w (Ising J_ij = -w)

struct Matrix {
    std::vector<int8_t> J = std::vector<int8_t>(size_t(sca::NP) * sca::NP, 0);   // J[x * NP + y]; symmetric, J_xx = 0, padding 0
    int64_t sumw = 0;          // sum_{i<j} w_ij with w = -J (cut = (sumw + sum_i s_i h_i / 2) / 2)
    int maxabs = 0;            // max |J|
    std::string source;        // description
    std::vector<Edge> edges;   // raw edges when loaded from an edge list (independent scorer); empty otherwise
    bool from_packed = false;  // K2000: the packed +-1 graph is kept for the original independent scorer
    int n = sca::N;            // r3: active spins (1..NP); rows/columns >= n are zero
    std::vector<int32_t> bias = std::vector<int32_t>(sca::NP, 0);   // r3: b_i (0 for padding)
    bool has_bias = false;     // r3: problem with a linear term (bias mode on the engine)
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
            if ((x == y || x >= m.n || y >= m.n) && v != 0) throw std::runtime_error("nonzero diagonal or padding");
            mx = std::max(mx, std::abs(v));
            if (x < y) s -= v;
        }
    m.sumw = s; m.maxabs = mx;
    for (int i = m.n; i < sca::NP; ++i) if (m.bias[i] != 0) throw std::runtime_error("nonzero padding bias");
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
    if (!(f >> n >> e) || n < 1 || n > sca::NP) throw std::runtime_error("G-set instance needs 1 <= N <= 2048: " + path);
    Matrix m; m.source = "gset:" + path; m.n = n;
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
inline Matrix random_k(uint64_t seed, double density, int maxabs, int n = sca::N, int bmax = 0) {
    Matrix m; m.n = n; m.source = "rand:" + std::to_string(seed) + ":" + std::to_string(density) + ":" + std::to_string(maxabs) +
                                   ":" + std::to_string(n) + ":" + std::to_string(bmax);
    uint64_t st = seed;
    const uint64_t thr = uint64_t(density * 18446744073709551616.0 - 1.0);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const uint64_t a = sca::Rng::splitmix(st), b = sca::Rng::splitmix(st);
            if (a > thr) continue;
            const int mag = 1 + int(b % uint64_t(maxabs)); const int v = ((b >> 32) & 1) ? mag : -mag;
            at(m, i, j) = int8_t(v); at(m, j, i) = int8_t(v);
            m.edges.push_back({i, j, -v});
        }
    if (bmax > 0) {   // r3: uniform integer biases in [-bmax, bmax]
        m.has_bias = true;
        for (int i = 0; i < n; ++i) m.bias[i] = int32_t(sca::Rng::splitmix(st) % uint64_t(2 * bmax + 1)) - bmax;
    }
    finish(m);
    return m;
}

// r3: text Ising problem. Lines: "n <n>", "J <i> <j> <value>" (0-based, i != j; symmetric), "b <i> <value>", '#' comments.
// H(s) = -sum_{i<j} J_ij s_i s_j - sum_i b_i s_i (MULTIBIT_SPEC.md). Also accepts "target <energy>" (recorded, not used here).
inline Matrix from_ising(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    Matrix m; m.source = "ising:" + path; m.n = 0; m.has_bias = true;
    std::string tok;
    while (f >> tok) {
        if (tok[0] == '#') { std::string rest; std::getline(f, rest); continue; }
        if (tok == "n") { f >> m.n; if (m.n < 1 || m.n > sca::NP) throw std::runtime_error("bad n"); }
        else if (tok == "J") { int i, j, v; f >> i >> j >> v;
            if (i < 0 || j < 0 || i >= m.n || j >= m.n || i == j || std::abs(v) > 127) throw std::runtime_error("bad J entry");
            at(m, i, j) = int8_t(v); at(m, j, i) = int8_t(v); if (i < j) m.edges.push_back({i, j, -v}); }
        else if (tok == "b") { int i; long v; f >> i >> v; if (i < 0 || i >= m.n) throw std::runtime_error("bad b entry"); m.bias[i] = int32_t(v); }
        else if (tok == "target") { std::string rest; std::getline(f, rest); }
        else throw std::runtime_error("bad token " + tok);
    }
    if (m.n == 0) throw std::runtime_error("missing n");
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
    for (int x = 0; x < m.n; ++x)
        for (int y = 0; y < m.n; ++y) { const double v = at(m, x, y); s += v * v; }
    return s / m.n / m.n;
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

// r3: field bound of the engine (MULTIBIT_SPEC widths): max_y sum_x |J_xy| + max_i |b_i| < 2^(HB-1); the host refuses otherwise.
inline int64_t field_bound(const Matrix& m) {
    int64_t mx = 0, bm = 0;
    for (int y = 0; y < sca::NP; ++y) { int64_t s = 0; for (int x = 0; x < sca::NP; ++x) s += std::abs(int(at(m, x, y))); mx = std::max(mx, s); }
    for (int i = 0; i < sca::NP; ++i) bm = std::max<int64_t>(bm, std::llabs(int64_t(m.bias[i])));
    return mx + bm;
}
// r3: bias stream segment, 256 beats of 128 bits (as 2 x 64-bit words lo, hi): beat i holds the biases of group i >> 3,
// fields j = 8 (i & 7) + q, q = 0..7 (16-bit two's complement lanes), i.e. spin y = (i >> 3) + 32 j.
inline std::vector<uint64_t> bias_beats(const Matrix& m) {
    std::vector<uint64_t> w(512, 0);
    for (int i = 0; i < 256; ++i)
        for (int q = 0; q < 8; ++q) {
            const int y = (i >> 3) + 32 * (8 * (i & 7) + q); const int32_t b = m.bias[y];
            if (b < -32768 || b > 32767) throw std::runtime_error("bias does not fit 16 bits");
            w[2 * i + (q >> 2)] |= uint64_t(uint16_t(int16_t(b))) << (16 * (q & 3));
        }
    return w;
}
// r3: independent scorer for problems with bias: sum_{i<n} s_i h_i with h_i = b_i + sum_j J_ij s_j, and sum_i b_i s_i.
inline std::pair<int64_t, int64_t> sh_bs_of(const Matrix& m, const uint64_t* bits) {
    auto sv = [&](int i) { return ((bits[i >> 6] >> (i & 63)) & 1) ? 1 : -1; };
    int64_t sh = 0, bs = 0;
    for (int i = 0; i < m.n; ++i) {
        int64_t h = m.bias[i];
        for (int j = 0; j < m.n; ++j) if (j != i) h += int64_t(at(m, i, j)) * sv(j);
        sh += sv(i) * h; bs += int64_t(m.bias[i]) * sv(i);
    }
    return {sh, bs};
}

inline Matrix load(const std::string& spec) {
    if (spec.rfind("k2000:", 0) == 0) return from_k2000(spec.substr(6));
    if (spec.rfind("gset:", 0) == 0) return from_gset(spec.substr(5));
    if (spec.rfind("rand:", 0) == 0) {   // rand:<seed>:<density>:<maxabs>[:<n>[:<bmax>]]
        std::stringstream ss(spec.substr(5)); std::string a, b, c, d, e;
        std::getline(ss, a, ':'); std::getline(ss, b, ':'); std::getline(ss, c, ':'); std::getline(ss, d, ':'); std::getline(ss, e, ':');
        return random_k(std::stoull(a), std::stod(b), std::stoi(c), d.empty() ? sca::N : std::stoi(d), e.empty() ? 0 : std::stoi(e));
    }
    if (spec.rfind("ising:", 0) == 0) return from_ising(spec.substr(6));
    throw std::runtime_error("matrix spec must be k2000:<file>, gset:<file>, ising:<file> or rand:<seed>:<density>:<maxabs>[:<n>[:<bmax>]]");
}

}  // namespace mb3
