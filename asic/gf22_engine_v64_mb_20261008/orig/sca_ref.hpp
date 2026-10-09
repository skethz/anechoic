// Bit-exact reference for the V80 Onsager-SCA engine (see ../DESIGN.md).
// The HLS kernel must reproduce these results exactly: final spins, cut, total flips, and the per-step n_lin trace.
#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace sca {

constexpr int N = 2000;      // active spins (K2000)
constexpr int NP = 2048;     // padded spins
constexpr int WORDS = NP / 64;
#ifndef SCA_LANES
#define SCA_LANES 128         // v4 engine (board image e818b91b...) uses 128; v5 builds pass -DSCA_LANES=256
#endif
constexpr int LANES = SCA_LANES;   // decision lanes; spin i is handled by lane i % LANES in round i / LANES
constexpr int ROUNDS = NP / LANES;

// Packed couplings. Bit (y % 64) of word J[x][y / 64] is 1 iff J_xy = +1. Diagonal and padding are masked separately.
struct Graph {
    std::vector<std::array<uint64_t, WORDS>> J = std::vector<std::array<uint64_t, WORDS>>(NP);
    int64_t sumw = 0;  // sum_{i<j} w_ij, with w = -J
};

inline int jval(const Graph& g, int x, int y) { return (g.J[x][y >> 6] >> (y & 63)) & 1 ? 1 : -1; }

inline Graph load_graph(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[8]; uint32_t n = 0;
    if (std::fread(magic, 1, 8, f) != 8 || std::string(magic, 8) != "SNOWGPU1" || std::fread(&n, 4, 1, f) != 1 || n != N)
        throw std::runtime_error("bad graph header");
    Graph g;
    for (int x = 0; x < N; ++x) {
        if (std::fread(g.J[x].data(), 8, WORDS, f) != WORDS) throw std::runtime_error("short graph file");
        g.J[x][WORDS - 1] &= (uint64_t(1) << (N - 64 * (WORDS - 1))) - 1;  // zero padding columns
        g.J[x][x >> 6] &= ~(uint64_t(1) << (x & 63));                      // stored diagonal bit is irrelevant; masked in use
    }
    std::fclose(f);
    int64_t s = 0;
    for (int x = 0; x < N; ++x)
        for (int y = x + 1; y < N; ++y) s -= jval(g, x, y);
    g.sumw = s;
    return g;
}

// Per-step host tables (fixed point). fourT: Q16.16 of 4T_t; q: Q16.16; kcorr: Q8.24 of lambda_t/(2 T_{t-1}); kconst: Q16.16 signed.
struct Tables {
    std::vector<int32_t> fourT, q, kconst;
    std::vector<int64_t> kcorr;
};

inline Tables make_tables(const std::vector<double>& T, const std::vector<double>& q, const std::vector<double>& lam,
                          double tec_jv = 0.0) {
    const size_t S = T.size();
    Tables tb; tb.fourT.resize(S); tb.q.resize(S); tb.kcorr.resize(S); tb.kconst.resize(S);
    for (size_t t = 0; t < S; ++t) {
        tb.fourT[t] = int32_t(std::llround(4.0 * T[t] * 65536.0));
        tb.q[t] = int32_t(std::llround(q[t] * 65536.0));
        tb.kcorr[t] = t == 0 ? 0 : int64_t(std::llround(lam[t] / (2.0 * T[t - 1]) * 16777216.0));
        tb.kconst[t] = t == 0 ? 0 : int32_t(std::llround(-tec_jv * 65536.0));  // TEC: field + Jv*s_prev == corr = -Jv; no s_prev at t=0
    }
    return tb;
}

// xoshiro128**, seeded by splitmix64(seed, trial, lane).
struct Rng {
    uint32_t s[4];
    static uint64_t splitmix(uint64_t& x) {
        uint64_t z = (x += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    void seed(uint64_t seed, uint32_t trial, uint32_t lane) {
        uint64_t x = seed ^ (uint64_t(trial) << 20) ^ (uint64_t(lane) << 48) ^ 0x5CA0F00Dull;
        uint64_t a = splitmix(x), b = splitmix(x);
        s[0] = uint32_t(a); s[1] = uint32_t(a >> 32); s[2] = uint32_t(b); s[3] = uint32_t(b >> 32);
        if (!(s[0] | s[1] | s[2] | s[3])) s[0] = 1;
    }
    static uint32_t rotl(uint32_t x, int k) { return (x << k) | (x >> (32 - k)); }
    uint32_t next() {
        uint32_t r = rotl(s[1] * 5u, 7) * 9u, t = s[1] << 9;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 11);
        return r;
    }
};

struct Result {
    std::array<uint64_t, WORDS> spins{};   // bit = 1 means +1
    int64_t cut = 0;
    int64_t flips = 0;
    std::vector<int32_t> n_lin;            // per step
};

// dense (optional): NP*NP int8 copy of J (0 on diagonal/padding) used only to speed up the same integer arithmetic.
inline Result run_trial(const Graph& g, const Tables& tb, uint64_t seed, uint32_t trial, bool trace = false,
                        const int8_t* dense = nullptr) {
    const int S = int(tb.fourT.size());
    std::vector<Rng> rng(LANES);
    for (int l = 0; l < LANES; ++l) rng[l].seed(seed, trial, l);
    std::vector<int8_t> s(NP, 1), sprev(NP, 1);
    for (int k = 0; k < ROUNDS; ++k)
        for (int l = 0; l < LANES; ++l) {
            int i = k * LANES + l; uint32_t r = rng[l].next();
            if (i < N) s[i] = (r >> 31) ? 1 : -1;
        }
    std::vector<int32_t> h(NP, 0);
    for (int x = 0; x < N; ++x)
        for (int y = 0; y < N; ++y)
            if (y != x) h[y] += (dense ? dense[size_t(x) * NP + y] : jval(g, x, y)) * s[x];
    Result res; if (trace) res.n_lin.resize(S);
    int32_t nlin_prev = 0; std::vector<int> flipped; flipped.reserve(NP);
    for (int t = 0; t < S; ++t) {
        const int64_t corr = ((int64_t(nlin_prev) * tb.kcorr[t]) >> 8) + tb.kconst[t];
        const int32_t fourT = tb.fourT[t], twoT = fourT >> 1;
        int32_t nlin = 0; flipped.clear();
        for (int k = 0; k < ROUNDS; ++k)
            for (int l = 0; l < LANES; ++l) {
                const int i = k * LANES + l;
                const uint32_t u16 = rng[l].next() >> 16;
                if (i >= N) continue;
                const int64_t z = (int64_t(s[i]) * h[i]) * 65536 + tb.q[t] - (s[i] == sprev[i] ? corr : -corr);
                const int64_t r = ((int64_t(u16) - 32768) * fourT) >> 16;   // (u - 1/2)*4T in [-2T, 2T); arithmetic shift (floor)
                if (z > -twoT && z < twoT) ++nlin;
                if (z < r) flipped.push_back(i);
            }
        for (int i = 0; i < N; ++i) sprev[i] = s[i];
        for (int x : flipped) s[x] = int8_t(-s[x]);
        for (int x : flipped) {
            const int d = 2 * s[x];
            if (dense) {
                const int8_t* row = dense + size_t(x) * NP;
                for (int y = 0; y < N; ++y) h[y] += d * row[y];   // row[x] == 0
            } else {
                for (int y = 0; y < N; ++y)
                    if (y != x) h[y] += d * jval(g, x, y);
            }
        }
        res.flips += int64_t(flipped.size());
        nlin_prev = nlin;
        if (trace) res.n_lin[t] = nlin;
    }
    int64_t sh = 0;
    for (int i = 0; i < N; ++i) {
        sh += int64_t(s[i]) * h[i];
        if (s[i] > 0) res.spins[i >> 6] |= uint64_t(1) << (i & 63);
    }
    res.cut = (g.sumw + sh / 2) / 2;
    return res;
}

inline std::vector<int8_t> dense_matrix(const Graph& g) {
    std::vector<int8_t> d(size_t(NP) * NP, 0);
    for (int x = 0; x < N; ++x)
        for (int y = 0; y < N; ++y)
            if (y != x) d[size_t(x) * NP + y] = int8_t(jval(g, x, y));
    return d;
}

}  // namespace sca
