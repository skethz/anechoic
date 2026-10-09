// Statistical check: the fixed-point engine at SCA_LANES (compile-time) on the O1 configuration, 2048 trials.
#include "sca_ref.hpp"
#include <atomic>
#include <cmath>
#include <thread>
using namespace sca;
int main(int argc, char** argv) {
    Graph g = load_graph(argv[1]); auto dense = dense_matrix(g);
    const int S = 960, TR = 2048; const double t0 = 12, q = 8, lam = 1.05;
    std::vector<double> T(S), Q(S, q), L(S, lam);
    for (int t = 0; t < S; ++t) { T[t] = t0 * std::pow(5.0 / t0, double(t) / (S - 1)); if (t >= 0.7 * S) L[t] = lam * (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S); }
    auto tb = make_tables(T, Q, L);
    std::atomic<int> next{0}, succ{0}; std::atomic<long long> fl{0}, cs{0};
    std::vector<std::thread> th;
    for (unsigned w = 0; w < std::thread::hardware_concurrency(); ++w) th.emplace_back([&] {
        for (int i; (i = next++) < TR;) { auto r = run_trial(g, tb, 20261005ull, uint32_t(i), false, dense.data()); succ += r.cut >= 33000; fl += r.flips; cs += r.cut; }
    });
    for (auto& x : th) x.join();
    double p = double(succ) / TR;
    std::printf("LANES=%d O1: p=%.3f (+-%.3f) mean_cut=%.1f mean_flips=%.0f\n", LANES, p, 1.96 * std::sqrt(p * (1 - p) / TR), double(cs) / TR, double(fl) / TR);
}
