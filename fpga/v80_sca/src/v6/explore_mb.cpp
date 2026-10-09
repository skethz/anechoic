// Reference-model schedule exploration for integer coupling matrices (sca_ref.hpp run_trial(..., dense), unchanged).
// Not part of any frozen protocol: used to choose test schedules for the RTL regression and to sanity-check G-set
// configurations. One configuration per invocation (parallelised by the caller).
// Build: c++ -O3 -std=c++17 -DSCA_LANES=256 explore_mb.cpp -o explore_mb
// Usage: explore_mb <matrix> <target> <trials> <seed> <offset> <t0> <t1> <S> <q> <lambda> <ramp 0/1> <jv> <kappa> <scale|auto> [name]
// Prints one JSON line: mean/max cut, success fraction, mean flips, v6.4 cycle-model cycles per trial.
#include "../mb_common.hpp"
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc < 15) { std::fprintf(stderr, "usage: see header\n"); return 2; }
    try {
        const mb::Matrix mat = mb::load(argv[1]);
        const long long target = std::atoll(argv[2]); const int trials = std::atoi(argv[3]);
        const uint64_t seed = std::strtoull(argv[4], 0, 0); const uint32_t off = uint32_t(std::strtoul(argv[5], 0, 0));
        mb::Sched c; c.t0 = std::atof(argv[6]); c.t1 = std::atof(argv[7]); c.S = std::atoi(argv[8]); c.q = std::atof(argv[9]);
        c.lam = std::atof(argv[10]); c.ramp = std::atoi(argv[11]) != 0; c.jv = std::atof(argv[12]); c.kappa = std::atof(argv[13]);
        c.scale = std::string(argv[14]) == "auto" ? mb::j2_scale(mat) : std::atof(argv[14]);
        const std::string name = argc > 15 ? argv[15] : "";
        const sca::Tables tb = mb::tables(c);
        mb::check_tables(tb, 4);
        const sca::Graph g = mb::graph_of(mat);
        double sc = 0, sf = 0; long long mx = -(1ll << 62); int hits = 0;
        for (int t = 0; t < trials; ++t) {
            const sca::Result r = sca::run_trial(g, tb, seed, off + uint32_t(t), false, mat.J.data());
            sc += double(r.cut); sf += double(r.flips); mx = std::max(mx, (long long)r.cut); hits += r.cut >= target;
        }
        const double flips = sf / trials, cyc = 0.1424 * flips + 18.09 * c.S + 886;
        std::printf("{\"name\":\"%s\",\"t0\":%g,\"t1\":%g,\"S\":%d,\"q\":%g,\"lambda\":%g,\"ramp\":%d,\"jv\":%g,\"kappa\":%g,\"scale\":%.9g,"
                    "\"trials\":%d,\"target\":%lld,\"hits\":%d,\"p\":%.6f,\"mean_cut\":%.3f,\"max_cut\":%lld,\"mean_flips\":%.1f,"
                    "\"model_cycles\":%.1f}\n", name.c_str(), c.t0, c.t1, c.S, c.q, c.lam, int(c.ramp), c.jv, c.kappa, c.scale, trials,
                    target, hits, double(hits) / trials, sc / trials, mx, flips, cyc);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL: %s\n", e.what());
        return 1;
    }
}
