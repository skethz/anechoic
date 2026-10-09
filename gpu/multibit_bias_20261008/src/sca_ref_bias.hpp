// Golden reference with a per-spin bias (external field) and a runtime number of active spins.
// Spec: ../MULTIBIT_SPEC.md, section "Bias (external field)". sca_ref.hpp is unchanged.
// This is sca::run_trial (dense path) with three differences:
//   (1) the initial fields include the bias: h_i(0) = b_i + sum_{j != i} J_ij s_j(0); the per-flip update is unchanged,
//       so the bias persists exactly;
//   (2) the number of active spins n (1..NP) is a runtime argument. Spins i >= n are padding: never updated (s = +1),
//       and the random-number streams are consumed exactly as in run_trial, independent of n;
//   (3) the couplings come only from the dense NP x NP int8 matrix (0 on the diagonal and padding), which carries
//       multi-bit couplings through the same integer arithmetic.
// With n = N, bias = nullptr (or all zeros) and the K2000 dense matrix, it returns what run_trial returns
// (see test_ref_bias.cpp).
#pragma once
#include "sca_ref.hpp"

namespace sca {

struct ResultBias {
    std::array<uint64_t, WORDS> spins{};   // bit = 1 means +1
    int64_t sum_sh = 0;                    // sum_i s_i h_i with the biased fields (what the device reports)
    int64_t sum_bs = 0;                    // sum_i b_i s_i, so that H = -(sum_sh + sum_bs) / 2
    int64_t flips = 0;
    std::vector<int32_t> n_lin;            // per step (filled when trace is set)
};

inline ResultBias run_trial_bias(int n, const int8_t* dense, const int32_t* bias, const Tables& tb, uint64_t seed,
                                 uint32_t trial, bool trace = false) {
    if (n < 1 || n > NP) throw std::runtime_error("run_trial_bias: n out of range");
    const int S = int(tb.fourT.size());
    std::vector<Rng> rng(LANES);
    for (int l = 0; l < LANES; ++l) rng[l].seed(seed, trial, l);
    std::vector<int8_t> s(NP, 1), sprev(NP, 1);
    for (int k = 0; k < ROUNDS; ++k)
        for (int l = 0; l < LANES; ++l) {
            const int i = k * LANES + l; const uint32_t r = rng[l].next();
            if (i < n) s[i] = (r >> 31) ? 1 : -1;
        }
    std::vector<int32_t> h(NP, 0);
    for (int y = 0; y < n; ++y) h[y] = bias ? bias[y] : 0;
    for (int x = 0; x < n; ++x)
        for (int y = 0; y < n; ++y)
            if (y != x) h[y] += int32_t(dense[size_t(x) * NP + y]) * s[x];
    ResultBias res; if (trace) res.n_lin.resize(S);
    int32_t nlin_prev = 0; std::vector<int> flipped; flipped.reserve(NP);
    for (int t = 0; t < S; ++t) {
        const int64_t corr = ((int64_t(nlin_prev) * tb.kcorr[t]) >> 8) + tb.kconst[t];
        const int32_t fourT = tb.fourT[t], twoT = fourT >> 1;
        int32_t nlin = 0; flipped.clear();
        for (int k = 0; k < ROUNDS; ++k)
            for (int l = 0; l < LANES; ++l) {
                const int i = k * LANES + l;
                const uint32_t u16 = rng[l].next() >> 16;
                if (i >= n) continue;
                const int64_t z = (int64_t(s[i]) * h[i]) * 65536 + tb.q[t] - (s[i] == sprev[i] ? corr : -corr);
                const int64_t r = ((int64_t(u16) - 32768) * fourT) >> 16;   // (u - 1/2)*4T in [-2T, 2T)
                if (z > -twoT && z < twoT) ++nlin;
                if (z < r) flipped.push_back(i);
            }
        for (int i = 0; i < n; ++i) sprev[i] = s[i];
        for (int x : flipped) s[x] = int8_t(-s[x]);
        for (int x : flipped) {
            const int d = 2 * s[x];
            const int8_t* row = dense + size_t(x) * NP;
            for (int y = 0; y < n; ++y) h[y] += d * row[y];   // row[x] == 0
        }
        res.flips += int64_t(flipped.size());
        nlin_prev = nlin;
        if (trace) res.n_lin[t] = nlin;
    }
    for (int i = 0; i < n; ++i) {
        res.sum_sh += int64_t(s[i]) * h[i];
        if (bias) res.sum_bs += int64_t(bias[i]) * s[i];
        if (s[i] > 0) res.spins[i >> 6] |= uint64_t(1) << (i & 63);
    }
    return res;
}

}  // namespace sca
