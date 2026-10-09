// Stream-level test vectors for the r3 core (sca_core_mb_r3.v): runtime active spins n and, in bias mode, per-spin biases.
// Expected outputs from the shared golden model sca_ref_bias.hpp run_trial_bias(n, dense, bias, tables, seed, trial, trace)
// (unchanged; with b = 0 and n = N it equals run_trial). Copy of gen_vectors_mb.cpp's stream logic with the r3 header
// fields: H0[47:36] = n, H0[35] = bias mode (then a 256-beat bias segment follows the couplings) and output word 32 =
// sum_i s_i h_i in bias mode, else the cut (sumw + sum_sh / 2) / 2 as before.
// Build: c++ -O2 -std=c++17 -DSCA_LANES=256 gen_vectors_mb3.cpp -o gen_vectors_mb3
// Usage: gen_vectors_mb3 <K> <matrix> <outdir> <cases file> [bias: 0 | 1 | auto (default: auto = matrix has a bias)] [n: 0 = write 0 (build-time N)]
//   matrix: k2000:<file> | gset:<file> | ising:<file> | rand:<seed>:<density>:<maxabs>[:<n>[:<bmax>]]
//   cases file lines: name t0 t1 S q lambda ramp jv kappa scale(1|auto) trials trace keep   ('#' comments)
#include "../mb_common3.hpp"
#include <cstdlib>
#include <fstream>

int main(int argc, char** argv) {
    using namespace sca;
    if (argc < 5) { std::fprintf(stderr, "usage: gen_vectors_mb3 K matrix outdir cases [bias 0|1|auto] [nfield]\n"); return 2; }
    try {
        const int K = std::atoi(argv[1]);
        const mb3::Matrix mat = mb3::load(argv[2]);
        const std::string dir = argv[3];
        const std::string bmode_arg = argc > 5 ? argv[5] : "auto";
        const bool bmode = bmode_arg == "auto" ? mat.has_bias : bmode_arg == "1";
        const int nfield = argc > 6 ? std::atoi(argv[6]) : mat.n;   // value written into H0[47:36] (0 = build-time N)
        if (nfield != 0 && nfield != mat.n) throw std::runtime_error("n field must be 0 or the matrix's n");
        if (nfield == 0 && mat.n != N) throw std::runtime_error("n field 0 means N = 2000");
        const std::vector<uint64_t> P = mb3::planes(mat, K);
        const std::vector<uint64_t> Bw = mb3::bias_beats(mat);
        const Graph g = mb3::graph_of(mat);
        const double auto_scale = mb3::j2_scale(mat);
        std::ifstream cf(argv[4]);
        if (!cf) throw std::runtime_error("cannot open cases file");
        struct Case { std::string name; mb3::Sched s; int trials, trace, keep; };
        std::vector<Case> cases; std::string line;
        while (std::getline(cf, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line); Case c; std::string scale; int ramp;
            ss >> c.name >> c.s.t0 >> c.s.t1 >> c.s.S >> c.s.q >> c.s.lam >> ramp >> c.s.jv >> c.s.kappa >> scale >> c.trials >> c.trace >> c.keep;
            if (!ss) throw std::runtime_error("bad cases line: " + line);
            c.s.ramp = ramp != 0; c.s.scale = scale == "auto" ? auto_scale : std::stod(scale);
            cases.push_back(c);
        }
        FILE* fs = std::fopen((dir + "/stim.hex").c_str(), "w");
        FILE* fe = std::fopen((dir + "/expect.hex").c_str(), "w");
        FILE* fc = std::fopen((dir + "/counts.txt").c_str(), "w");
        FILE* fi = std::fopen((dir + "/cases.txt").c_str(), "w");
        if (!fs || !fe || !fc || !fi) { std::perror("open"); return 1; }
        auto beat = [&](uint64_t hi, uint64_t lo) { std::fprintf(fs, "%016llx%016llx\n", (unsigned long long)hi, (unsigned long long)lo); };
        auto word = [&](uint64_t v) { std::fprintf(fe, "%016llx\n", (unsigned long long)v); };
        std::fprintf(fi, "K %d matrix %s n %d nfield %d bias_mode %d maxabs %d sumw %lld field_bound %lld fingerprint %016llx auto_scale %.12g\n",
                     K, mat.source.c_str(), mat.n, nfield, int(bmode), mat.maxabs, (long long)mat.sumw, (long long)mb3::field_bound(mat),
                     (unsigned long long)mb3::fingerprint(mat), auto_scale);
        long nbeats = 0, nwords = 0;
        for (size_t ci = 0; ci < cases.size(); ++ci) {
            const Case& c = cases[ci];
            const Tables tb = mb3::tables(c.s);
            mb3::check_tables(tb, K);
            const uint64_t seed = 0xC0FFEEull + ci; const uint32_t off = 1000 + 100 * uint32_t(ci);
            const int loadJ = ci == 0;
            beat(seed, (uint64_t(nfield & 0xFFF) << 36) | (uint64_t(bmode) << 35) | (uint64_t(c.keep) << 34) | (uint64_t(c.trace) << 33) |
                       (uint64_t(loadJ) << 32) | (uint64_t(c.s.S) << 16) | uint64_t(c.trials));
            beat(uint64_t(g.sumw), uint64_t(off));
            nbeats += 2;
            if (loadJ) {
                for (size_t r = 0; r < size_t(K) * NP; ++r)
                    for (int cc = 0; cc < 16; ++cc) beat(P[r * WORDS + 2 * cc + 1], P[r * WORDS + 2 * cc]);
                nbeats += long(K) * NP * 16;
                if (bmode) { for (int i = 0; i < 256; ++i) beat(Bw[2 * i + 1], Bw[2 * i]); nbeats += 256; }
            }
            for (int t = 0; t < (c.keep ? 0 : c.s.S); ++t)
                beat((uint64_t(uint32_t(tb.kcorr[t])) << 32) | uint64_t(uint32_t(tb.kconst[t])),
                     (uint64_t(uint32_t(tb.q[t])) << 32) | uint64_t(uint32_t(tb.fourT[t])));
            nbeats += c.keep ? 0 : c.s.S;
            std::fprintf(fi, "case %zu %s S %d t0 %g t1 %g q %g lambda %g ramp %d jv %g kappa %g scale %.12g trials %d trace %d keep %d seed %llu offset %u\n",
                         ci, c.name.c_str(), c.s.S, c.s.t0, c.s.t1, c.s.q, c.s.lam, int(c.s.ramp), c.s.jv, c.s.kappa, c.s.scale, c.trials,
                         c.trace, c.keep, (unsigned long long)seed, off);
            for (int tr = 0; tr < c.trials; ++tr) {
                const ResultBias r = run_trial_bias(mat.n, mat.J.data(), bmode ? mat.bias.data() : nullptr, tb, seed, off + tr, true);
                const int64_t cut = (g.sumw + r.sum_sh / 2) / 2;
                for (int w = 0; w < WORDS; ++w) word(r.spins[w]);
                word(uint64_t(bmode ? r.sum_sh : cut)); word(uint64_t(r.flips));
                nwords += WORDS + 2;
                if (c.trace) { for (int t = 0; t < c.s.S; ++t) word(uint64_t(uint32_t(r.n_lin[t]))); nwords += c.s.S; }
                // independent checks of the golden result (scorers from the raw problem)
                if (bmode) {
                    const auto shbs = mb3::sh_bs_of(mat, r.spins.data());
                    if (shbs.first != r.sum_sh || shbs.second != r.sum_bs) throw std::runtime_error("independent sum_sh/sum_bs disagree");
                } else if (mb3::cut_of(mat, r.spins.data()) != cut) throw std::runtime_error("independent cut disagrees");
                std::printf("%-22s trial %d: %s %lld flips %lld energy %lld\n", c.name.c_str(), tr, bmode ? "sum_sh" : "cut",
                            (long long)(bmode ? r.sum_sh : cut), (long long)r.flips, (long long)(-(r.sum_sh + r.sum_bs) / 2));
                std::fprintf(fi, "  trial %u %s %lld flips %lld energy2 %lld\n", off + tr, bmode ? "sum_sh" : "cut",
                             (long long)(bmode ? r.sum_sh : cut), (long long)r.flips, (long long)(-(r.sum_sh + r.sum_bs)));
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
