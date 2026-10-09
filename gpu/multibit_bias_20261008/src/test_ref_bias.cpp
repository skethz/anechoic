// Checks for sca_ref_bias.hpp.
// (1) b = 0, n = N on K2000: run_trial_bias equals run_trial bit for bit (spins, cut, flips, n_lin trace), for an
//     Onsager-online and a TEC-style schedule.
// (2) Random bias, random K-bit couplings, n < N: the final fields equal a from-scratch recomputation
//     h_i = b_i + sum_j J_ij s_j, so the bias persists through the incremental updates.
// Build: c++ -O2 -std=c++17 -DSCA_LANES=256 test_ref_bias.cpp -o test_ref_bias ; run: ./test_ref_bias <K2000.bin>
#include "sca_ref_bias.hpp"
#include <cstdio>
#include <random>

using namespace sca;

static std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S);
    for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, S > 1 ? double(t) / (S - 1) : 0.0);
    return T;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s K2000.bin\n", argv[0]); return 2; }
    const Graph g = load_graph(argv[1]);
    std::vector<int8_t> dense(size_t(NP) * NP, 0);
    for (int x = 0; x < N; ++x)
        for (int y = 0; y < N; ++y)
            if (x != y) dense[size_t(x) * NP + y] = int8_t(jval(g, x, y));
    int fails = 0, checks = 0;
    struct Sched { const char* name; double t0, t1, q, lam, jv; int S; };
    const Sched sch[] = {{"onsager-online", 12, 5, 6, 0.9, 0, 360}, {"tec", 30, 5, 8, 0, -4, 300}};
    for (const auto& sc : sch) {
        const int S = sc.S;
        const Tables tb = make_tables(geom(sc.t0, sc.t1, S), std::vector<double>(S, sc.q), std::vector<double>(S, sc.lam), sc.jv);
        for (uint32_t trial = 0; trial < 6; ++trial) {
            const Result a = run_trial(g, tb, 20261008, trial, true, dense.data());
            const ResultBias b = run_trial_bias(N, dense.data(), nullptr, tb, 20261008, trial, true);
            const int64_t cut_b = (g.sumw + b.sum_sh / 2) / 2;
            const bool ok = a.spins == b.spins && a.flips == b.flips && a.n_lin == b.n_lin && a.cut == cut_b;
            ++checks; fails += !ok;
            std::printf("[b=0, K2000] %-15s trial %u: cut %lld/%lld flips %lld/%lld -> %s\n", sc.name, trial,
                        (long long)a.cut, (long long)cut_b, (long long)a.flips, (long long)b.flips, ok ? "identical" : "MISMATCH");
        }
    }
    // (2) random K = 4 couplings (|J| <= 7), random bias, n = 800
    std::mt19937_64 gen(7);
    const int n = 800;
    std::vector<int8_t> dj(size_t(NP) * NP, 0);
    std::vector<int32_t> bias(NP, 0);
    std::uniform_int_distribution<int> uj(-7, 7), ub(-300, 300);
    for (int x = 0; x < n; ++x)
        for (int y = x + 1; y < n; ++y) { const int v = uj(gen); dj[size_t(x) * NP + y] = dj[size_t(y) * NP + x] = int8_t(v); }
    for (int i = 0; i < n; ++i) bias[i] = ub(gen);
    const int S = 400;
    const Tables tb = make_tables(geom(200, 20, S), std::vector<double>(S, 40), std::vector<double>(S, 0.0));
    for (uint32_t trial = 0; trial < 4; ++trial) {
        const ResultBias r = run_trial_bias(n, dj.data(), bias.data(), tb, 99, trial, true);
        std::vector<int> s(n);
        for (int i = 0; i < n; ++i) s[i] = (r.spins[i >> 6] >> (i & 63)) & 1 ? 1 : -1;
        int64_t sum_sh = 0, sum_bs = 0;
        for (int i = 0; i < n; ++i) {
            int64_t hi = bias[i];
            for (int j = 0; j < n; ++j) if (j != i) hi += int64_t(dj[size_t(i) * NP + j]) * s[j];
            sum_sh += s[i] * hi; sum_bs += int64_t(bias[i]) * s[i];
        }
        const bool ok = sum_sh == r.sum_sh && sum_bs == r.sum_bs;
        ++checks; fails += !ok;
        std::printf("[bias, K=4, n=%d] trial %u: sum_sh %lld/%lld, flips %lld -> %s\n", n, trial, (long long)r.sum_sh,
                    (long long)sum_sh, (long long)r.flips, ok ? "consistent" : "MISMATCH");
    }
    std::printf("%d of %d checks passed\n", checks - fails, checks);
    return fails ? 1 : 0;
}
