// M1 quality check: does the fixed-point/RNG engine reproduce the float model's success? (and dense==packed equivalence)
#include "sca_ref.hpp"
#include <atomic>
#include <cmath>
#include <thread>
using namespace sca;

static std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S); for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1)); return T;
}
static std::vector<double> lam_sched(double lam, bool ramp, int S) {
    std::vector<double> L(S, lam);
    if (ramp) for (int t = 0; t < S; ++t) if (t >= 0.7 * S) L[t] = lam * (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S);
    return L;
}
int main(int argc, char** argv) {
    Graph g = load_graph(argc > 1 ? argv[1] : "../../v80_snowball/data/K2000.bin");
    std::printf("sumw=%lld (expect -1040)\n", (long long)g.sumw);
    auto dense = dense_matrix(g);
    struct Cfg { const char* name; double t0; int S; double q; double lam; bool ramp; double py_p; };
    std::vector<Cfg> cfgs = {
        {"plain_T40_S1560_q4", 40, 1560, 4, 0.0, false, 0.824},      // Python held-out (A2): 0.824
        {"onsager0.75_T20_S960_q4", 20, 960, 4, 0.75, false, 0.866},  // Python held-out (A, E=1): 0.866
        {"onsager0.7ramp_T30_S560_q4", 30, 560, 4, 0.7, true, 0.699}, // Python held-out (A2): 0.699
    };
    // Equivalence: packed vs dense path, one short trace
    {
        auto tb = make_tables(geom(20, 5, 40), std::vector<double>(40, 4.0), lam_sched(0.7, true, 40));
        auto a = run_trial(g, tb, 77, 3, true), b = run_trial(g, tb, 77, 3, true, dense.data());
        std::printf("packed==dense: %s (cut %lld, flips %lld)\n", (a.spins == b.spins && a.n_lin == b.n_lin && a.cut == b.cut) ? "PASS" : "FAIL",
                    (long long)a.cut, (long long)a.flips);
    }
    const int TR = argc > 2 ? std::atoi(argv[2]) : 1024;
    for (auto& c : cfgs) {
        auto tb = make_tables(geom(c.t0, 5, c.S), std::vector<double>(c.S, c.q), lam_sched(c.lam, c.ramp, c.S));
        std::atomic<int> next{0}, succ{0}; std::atomic<long long> cutsum{0}, flipsum{0};
        unsigned nt = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> th;
        for (unsigned w = 0; w < nt; ++w) th.emplace_back([&] {
            for (int i; (i = next++) < TR;) {
                auto r = run_trial(g, tb, 20261004ull, uint32_t(i), false, dense.data());
                if (r.cut >= 33000) ++succ; cutsum += r.cut; flipsum += r.flips;
            }
        });
        for (auto& x : th) x.join();
        double p = double(succ) / TR, se = std::sqrt(p * (1 - p) / TR);
        std::printf("%-28s p=%.3f (+-%.3f)  python=%.3f  mean_cut=%.1f  mean_flips=%.0f\n", c.name, p, 1.96 * se, c.py_p,
                    double(cutsum) / TR, double(flipsum) / TR);
    }
}
