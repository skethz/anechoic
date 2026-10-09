// Stream-level test vectors for the multi-bit core (sca_core_mb.v): input beats and the expected output words, from the
// unchanged golden model sca_ref.hpp run_trial(..., dense). Copy of gen_vectors.cpp with K coupling planes (plane-major,
// see sca_core_mb.v) and integer matrices (../mb_common.hpp).
// Build: c++ -O2 -std=c++17 -DSCA_LANES=256 gen_vectors_mb.cpp -o gen_vectors_mb
// Usage: gen_vectors_mb <K> <matrix> <outdir> [suite [case]]   (K = 1: v6.4b packed format, K2000 only)
//   matrix: k2000:<K2000.bin> | gset:<file> | rand:<seed>:<density>:<maxabs>
//   suite : short (default: v6 short suite incl. resident tables; for k2000 the exact v6.4 cases) | proto (k2000 full
//           O4, X5, O1, P1, T1) | long (one full-length Onsager and TEC-T schedule for non-K2000 matrices)
//   case  : run only this case index of the suite (couplings loaded)
#include "../mb_common.hpp"
#include <cstdlib>

int main(int argc, char** argv) {
    using namespace sca;
    if (argc < 4) { std::fprintf(stderr, "usage: gen_vectors_mb K matrix outdir [suite [case]]\n"); return 2; }
    try {
        const int K = std::atoi(argv[1]);
        const mb::Matrix mat = mb::load(argv[2]);
        const std::string dir = argv[3];
        const std::string suite = argc > 4 ? argv[4] : "short";
        const bool k2000 = mat.from_packed;
        const double auto_scale = mb::j2_scale(mat);
        long nnz = 0; for (int8_t v : mat.J) nnz += v != 0;
        const double qd = std::round(2.0 * (0.3 * double(nnz) / N + 0.8)) / 2.0;   // ternary test suites only
        // K = 1: the v6.4b format (one packed +-1 plane, bit = 1 iff J = +1, stored diagonal bit 0), K2000 only, so that
        // v6.4b can be simulated on exactly the same cases (cycle comparison)
        if (K == 1 && !k2000) throw std::runtime_error("K = 1 (v6.4b format) needs the K2000 matrix");
        std::vector<uint64_t> P;
        if (K == 1) { for (int x = 0; x < NP; ++x) for (int w = 0; w < WORDS; ++w) P.push_back(mat.packed.J[x][w]); }
        else P = mb::planes(mat, K);
        const Graph g = mb::graph_of(mat);
        FILE* fs = std::fopen((dir + "/stim.hex").c_str(), "w");
        FILE* fe = std::fopen((dir + "/expect.hex").c_str(), "w");
        FILE* fc = std::fopen((dir + "/counts.txt").c_str(), "w");
        FILE* fi = std::fopen((dir + "/cases.txt").c_str(), "w");
        if (!fs || !fe || !fc || !fi) { std::perror("open"); return 1; }
        auto beat = [&](uint64_t hi, uint64_t lo) { std::fprintf(fs, "%016llx%016llx\n", (unsigned long long)hi, (unsigned long long)lo); };
        auto word = [&](uint64_t v) { std::fprintf(fe, "%016llx\n", (unsigned long long)v); };

        struct Case { const char* name; mb::Sched s; int trials; int trace; int keep; };
        auto S = [](double t0, double t1, int steps, double q, double lam, bool ramp, double jv, double kappa, double scale) {
            mb::Sched c; c.t0 = t0; c.t1 = t1; c.S = steps; c.q = q; c.lam = lam; c.ramp = ramp; c.jv = jv; c.kappa = kappa;
            c.scale = scale; return c;
        };
        std::vector<Case> cases;
        if (k2000 && suite == "short")   // exactly the v6.4 suite (gen_vectors.cpp): same tables, seeds and trial ids
            cases = {{"plain_S40", S(30, 5, 40, 4, 0, false, 0, 0, 1), 3, 1, 0},
                     {"onsager_ramp_S40", S(20, 5, 40, 4, 0.7, true, 0, 0, 1), 3, 1, 0},
                     {"tec_S40", S(20, 5, 40, 8, 0, false, -8, 0, 1), 2, 0, 0},
                     {"tec_S40_keep", S(20, 5, 40, 8, 0, false, -8, 0, 1), 2, 1, 1}};
        else if (k2000 && suite == "proto")   // gen_vectors.cpp "proto" (same order, seeds, trial ids), plus X5 as case 6
            cases = {{"O1", S(12, 5, 960, 8, 1.05, true, 0, 0, 1), 1, 0, 0},
                     {"O2", S(15, 5, 560, 6, 0.9, true, 0, 0, 1), 1, 0, 0},
                     {"O3", S(12, 5, 560, 6, 0.9, true, 0, 0, 1), 1, 0, 0},
                     {"O4", S(12, 5, 360, 6, 0.9, false, 0, 0, 1), 1, 0, 0},
                     {"P1", S(30, 5, 1560, 4, 0, false, 0, 0, 1), 1, 0, 0},
                     {"T1", S(30, 5, 1560, 8, 0, false, -4, 0, 1), 1, 0, 0},
                     {"X5", S(15, 5, 280, 8, 0, true, 0, 1.75, 1), 1, 0, 0}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "short")   // ternary (G-set): fields are O(degree), not O(N); the
            // self-coupling scales with the mean degree d (q = 0.3 d + 0.8, rounded to 0.5: G22 7, G32 2), which avoids the
            // all-equal period-2 attractor that small q gives on +1-weight graphs (explore_mb, build/explore_mb_20261007a)
            cases = {{"plain_S40", S(4, 0.5, 40, qd, 0, false, 0, 0, 1), 3, 1, 0},
                     {"onsager_ramp_S40", S(4, 0.5, 40, qd, 0.9, true, 0, 0, auto_scale), 3, 1, 0},
                     {"tec_S40", S(4, 0.5, 40, qd, 0, false, -1, 0, 1), 2, 0, 0},
                     {"tec_S40_keep", S(4, 0.5, 40, qd, 0, false, -1, 0, 1), 2, 1, 1}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "neg")   // r2 regression: small q drives +1-weight graphs into the
            // all-equal period-2 states (cut 0, sum_i s_i h_i < 0), the case whose device cut word was wrong in sca_core_mb.v
            cases = {{"plain_q1_S40", S(3, 0.3, 40, 1, 0, false, 0, 0, 1), 3, 1, 0},
                     {"onsager_q1_S40", S(2, 0.3, 40, 1, 0.9, true, 0, 0, auto_scale), 3, 1, 0},
                     {"tec_q1_S40", S(2, 0.3, 40, 1, 0, false, -0.5, 0, 1), 2, 0, 0}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "long")
            cases = {{"onsager_ramp_S600", S(4, 0.5, 600, qd, 0.9, true, 0, 0, auto_scale), 2, 1, 0},
                     {"tecT_ramp_S600", S(4, 0.5, 600, qd, 0, true, 0, 1.5, auto_scale), 1, 1, 0}};
        else if (!k2000 && suite == "short")   // multi-bit random matrix (e.g. K = 4, |J| <= 7)
            cases = {{"plain_S40", S(150, 20, 40, 20, 0, false, 0, 0, 1), 3, 1, 0},
                     {"onsager_ramp_S40", S(100, 20, 40, 20, 0.7, true, 0, 0, auto_scale), 3, 1, 0},
                     {"tec_S40", S(100, 20, 40, 20, 0, false, -20, 0, 1), 2, 0, 0},
                     {"tec_S40_keep", S(100, 20, 40, 20, 0, false, -20, 0, 1), 2, 1, 1}};
        else if (!k2000 && suite == "long")
            cases = {{"onsager_ramp_S400", S(150, 10, 400, 10, 0.9, true, 0, 0, auto_scale), 2, 1, 0},
                     {"tecT_ramp_S300", S(150, 10, 300, 20, 0, true, 0, 1.5, auto_scale), 1, 1, 0}};
        else throw std::runtime_error("unknown suite " + suite);
        if (argc > 5) { Case one = cases.at(std::atoi(argv[5])); one.keep = 0; cases = {one}; }

        std::fprintf(fi, "K %d matrix %s maxabs %d sumw %lld fingerprint %016llx auto_scale %.12g\n", K, mat.source.c_str(),
                     mat.maxabs, (long long)mat.sumw, (unsigned long long)mb::fingerprint(mat), auto_scale);
        long nbeats = 0, nwords = 0;
        for (size_t ci = 0; ci < cases.size(); ++ci) {
            const Case& c = cases[ci];
            const Tables tb = mb::tables(c.s);
            mb::check_tables(tb, K < 2 ? 2 : K);
            const uint64_t seed = 0xC0FFEEull + ci; const uint32_t off = 1000 + 100 * uint32_t(ci);
            const int loadJ = ci == 0;
            beat(seed, (uint64_t(c.keep) << 34) | (uint64_t(c.trace) << 33) | (uint64_t(loadJ) << 32) | (uint64_t(c.s.S) << 16) | uint64_t(c.trials));
            beat(uint64_t(g.sumw), uint64_t(off));
            nbeats += 2;
            if (loadJ) {
                for (size_t r = 0; r < size_t(K) * NP; ++r)
                    for (int cc = 0; cc < 16; ++cc) beat(P[r * WORDS + 2 * cc + 1], P[r * WORDS + 2 * cc]);
                nbeats += long(K) * NP * 16;
            }
            for (int t = 0; t < (c.keep ? 0 : c.s.S); ++t)
                beat((uint64_t(uint32_t(tb.kcorr[t])) << 32) | uint64_t(uint32_t(tb.kconst[t])),
                     (uint64_t(uint32_t(tb.q[t])) << 32) | uint64_t(uint32_t(tb.fourT[t])));
            nbeats += c.keep ? 0 : c.s.S;
            std::fprintf(fi, "case %zu %s S %d t0 %g t1 %g q %g lambda %g ramp %d jv %g kappa %g scale %.12g trials %d trace %d keep %d seed %llu offset %u\n",
                         ci, c.name, c.s.S, c.s.t0, c.s.t1, c.s.q, c.s.lam, int(c.s.ramp), c.s.jv, c.s.kappa, c.s.scale, c.trials,
                         c.trace, c.keep, (unsigned long long)seed, off);
            for (int tr = 0; tr < c.trials; ++tr) {
                const Result r = run_trial(g, tb, seed, off + tr, true, mat.J.data());
                for (int w = 0; w < WORDS; ++w) word(r.spins[w]);
                word(uint64_t(r.cut)); word(uint64_t(r.flips));
                nwords += WORDS + 2;
                if (c.trace) { for (int t = 0; t < c.s.S; ++t) word(uint64_t(uint32_t(r.n_lin[t]))); nwords += c.s.S; }
                const int64_t ind = mb::cut_of(mat, r.spins.data());
                if (ind != r.cut) throw std::runtime_error("independent cut disagrees with the reference");
                std::printf("%-20s trial %d: cut %lld flips %lld\n", c.name, tr, (long long)r.cut, (long long)r.flips);
                std::fprintf(fi, "  trial %u cut %lld flips %lld\n", off + tr, (long long)r.cut, (long long)r.flips);
            }
        }
        std::fprintf(fc, "%ld %ld\n", nbeats, nwords);
        std::fclose(fs); std::fclose(fe); std::fclose(fc); std::fclose(fi);
        std::printf("beats %ld words %ld\n", nbeats, nwords);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL: %s\n", e.what());
        return 1;
    }
}
