// Exploratory (not part of PROTOCOL_HW): does the bit-exact fixed-point engine reproduce the float model's p for the
// hardware cohort's largest deviation (P1), on fresh seeds? P2 and T1 are controls. Same reference the board matches bit for bit.
#include "sca_ref.hpp"
#include <atomic>
#include <cmath>
#include <thread>
using namespace sca;

static std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S); for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1)); return T;
}
int main(int argc, char** argv) {
    Graph g = load_graph(argv[1]);
    const uint64_t seed = std::strtoull(argv[2], nullptr, 10);
    const int TR = std::atoi(argv[3]);
    auto dense = dense_matrix(g);
    struct Cfg { const char* name; double t0; int S; double q; double jv; double p_J, p_hw; };
    std::vector<Cfg> cfgs = {{"P1 plain q4 T30 S1560", 30, 1560, 4, 0, 0.675, 0.613},
                             {"P2 plain q6 T30 S1560", 30, 1560, 6, 0, 0.567, 0.580},
                             {"T1 tec q8 Jv-4 T30 S1560", 30, 1560, 8, -4, 0.675, 0.631}};
    for (auto& c : cfgs) {
        auto tb = make_tables(geom(c.t0, 5, c.S), std::vector<double>(c.S, c.q), std::vector<double>(c.S, 0.0), c.jv);
        std::atomic<int> next{0}, succ{0}; std::atomic<long long> cutsum{0};
        std::vector<std::thread> th;
        for (unsigned w = 0; w < std::max(1u, std::thread::hardware_concurrency()); ++w) th.emplace_back([&] {
            for (int i; (i = next++) < TR;) {
                auto r = run_trial(g, tb, seed, uint32_t(i), false, dense.data());
                if (r.cut >= 33000) ++succ; cutsum += r.cut;
            }
        });
        for (auto& x : th) x.join();
        double p = double(succ) / TR, se = std::sqrt(p * (1 - p) / TR);
        std::printf("%-26s seed %llu  fixed-point p=%.3f (+-%.3f, n=%d)  float p_J=%.3f  board p_hw=%.3f  mean_cut=%.1f\n", c.name,
                    (unsigned long long)seed, p, 1.96 * se, TR, c.p_J, c.p_hw, double(cutsum) / TR);
        std::fflush(stdout);
    }
}
