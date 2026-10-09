// Stream-level test vectors for the v6 core (sca_core.v): input beats and the expected output words, from sca_ref.hpp.
// Build: c++ -O2 -std=c++17 -DSCA_LANES=256 gen_vectors.cpp -o gen_vectors
// Usage: gen_vectors <K2000.bin> <outdir> [full]
#include "../sca_ref.hpp"
#include <cmath>
#include <cstdlib>

static std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S);
    for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1));
    return T;
}

int main(int argc, char** argv) {
    using namespace sca;
    if (argc < 3) { std::fprintf(stderr, "usage: gen_vectors K2000.bin outdir [full]\n"); return 2; }
    const std::string dir = argv[2];
    const bool full = argc > 3;
    Graph g = load_graph(argv[1]);
    std::vector<int8_t> dense = dense_matrix(g);
    FILE* fs = std::fopen((dir + "/stim.hex").c_str(), "w");
    FILE* fe = std::fopen((dir + "/expect.hex").c_str(), "w");
    FILE* fc = std::fopen((dir + "/counts.txt").c_str(), "w");
    if (!fs || !fe || !fc) { std::perror("open"); return 1; }
    auto beat = [&](uint64_t hi, uint64_t lo) { std::fprintf(fs, "%016llx%016llx\n", (unsigned long long)hi, (unsigned long long)lo); };
    auto word = [&](uint64_t v) { std::fprintf(fe, "%016llx\n", (unsigned long long)v); };

    struct Case { const char* name; double t0; int S; double q; double lam; bool ramp; double jv; int trials; int trace; int keep = 0; };
    std::vector<Case> cases = {{"plain_S40", 30, 40, 4, 0, false, 0, 3, 1},
                               {"onsager_ramp_S40", 20, 40, 4, 0.7, true, 0, 3, 1},
                               {"tec_S40", 20, 40, 8, 0, false, -8, 2, 0},
                               {"tec_S40_keep", 20, 40, 8, 0, false, -8, 2, 1, 1}};   // v6.2: same tables, kept resident
    if (full) cases.push_back({"onsager_ramp_S560", 30, 560, 4, 0.7, true, 0, 2, 1});
    if (full && std::string(argv[3]) == "proto")   // frozen PROTOCOL_HW_MULTI configurations, one trial each, no trace
        cases = {{"O1", 12, 960, 8, 1.05, true, 0, 1, 0}, {"O2", 15, 560, 6, 0.9, true, 0, 1, 0},
                 {"O3", 12, 560, 6, 0.9, true, 0, 1, 0},  {"O4", 12, 360, 6, 0.9, false, 0, 1, 0},
                 {"P1", 30, 1560, 4, 0, false, 0, 1, 0},   {"T1", 30, 1560, 8, 0, false, -4, 1, 0}};
    if (argc > 4) { const Case one = cases.at(std::atoi(argv[4])); cases = {one}; }   // single case (couplings loaded)
    long nbeats = 0, nwords = 0;
    for (size_t ci = 0; ci < cases.size(); ++ci) {
        const Case& c = cases[ci];
        std::vector<double> T = geom(c.t0, 5, c.S), q(c.S, c.q), L(c.S, c.lam);
        if (c.ramp)
            for (int t = 0; t < c.S; ++t)
                if (t >= 0.7 * c.S) L[t] = c.lam * (c.S - 1 - t) / std::max(1.0, c.S - 1 - 0.7 * c.S);
        Tables tb = make_tables(T, q, L, c.jv);
        const uint64_t seed = 0xC0FFEEull + ci; const uint32_t off = 1000 + 100 * uint32_t(ci);
        const int loadJ = ci == 0;
        beat(seed, (uint64_t(c.keep) << 34) | (uint64_t(c.trace) << 33) | (uint64_t(loadJ) << 32) | (uint64_t(c.S) << 16) | uint64_t(c.trials));
        beat(uint64_t(g.sumw), uint64_t(off));
        nbeats += 2;
        if (loadJ) {
            for (int x = 0; x < NP; ++x)
                for (int cc = 0; cc < 16; ++cc) beat(g.J[x][2 * cc + 1], g.J[x][2 * cc]);
            nbeats += NP * 16;
        }
        for (int t = 0; t < (c.keep ? 0 : c.S); ++t) {
            if (tb.kcorr[t] != int64_t(int32_t(tb.kcorr[t]))) { std::fprintf(stderr, "kcorr does not fit int32\n"); return 1; }
            beat((uint64_t(uint32_t(tb.kcorr[t])) << 32) | uint64_t(uint32_t(tb.kconst[t])),
                 (uint64_t(uint32_t(tb.q[t])) << 32) | uint64_t(uint32_t(tb.fourT[t])));
        }
        nbeats += c.keep ? 0 : c.S;
        for (int tr = 0; tr < c.trials; ++tr) {
            Result r = run_trial(g, tb, seed, off + tr, true, dense.data());
            for (int w = 0; w < WORDS; ++w) word(r.spins[w]);
            word(uint64_t(r.cut)); word(uint64_t(r.flips));
            nwords += WORDS + 2;
            if (c.trace) { for (int t = 0; t < c.S; ++t) word(uint64_t(uint32_t(r.n_lin[t]))); nwords += c.S; }
            std::printf("%-20s trial %d: cut %lld flips %lld\n", c.name, tr, (long long)r.cut, (long long)r.flips);
        }
    }
    std::fprintf(fc, "%ld %ld\n", nbeats, nwords);
    std::fclose(fs); std::fclose(fe); std::fclose(fc);
    std::printf("beats %ld words %ld\n", nbeats, nwords);
    return 0;
}
