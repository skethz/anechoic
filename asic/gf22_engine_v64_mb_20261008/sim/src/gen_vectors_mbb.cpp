// Stream-level test vectors for the multi-bit ASIC core with bias and runtime n (rtl/sca_core_mb_asic.v).
// Copy of the FPGA gen_vectors_mb.cpp (orig/, unchanged suites and expected outputs from sca_ref.hpp run_trial(..., dense)
// for b = 0) plus a "bias" suite whose expected outputs come from sca_ref_bias.hpp run_trial_bias(n, dense, bias, ...).
// Build: c++ -O2 -std=c++17 -DSCA_LANES=256 gen_vectors_mbb.cpp -o gen_vectors_mbb
// Usage: gen_vectors_mbb <K> <matrix> <outdir> [suite [case]]
//   matrix: k2000:<K2000.bin> | gset:<file> | rand:<seed>:<density>:<maxabs> | randn:<seed>:<density>:<maxabs>:<n>
//           (randn: random matrix on the first n spins only, padding rows/columns zero)
//   suite : short | proto | long | neg (as gen_vectors_mb.cpp) | bias (randn matrices: random biases, runtime n,
//           raw score word; one case without bias)
// Stream additions of the ASIC core: H0 bit 35 = bias segment (256 beats after the coupling planes; beat a = int16
// biases 8a..8a+7, bias 8a+i in bits 16i..16i+15), bit 36 = biased field initialisation, bit 37 = raw score word
// (sum_i s_i h_i instead of the cut); H1[43:32] = active spins n (0 = N).
#include "mb_common.hpp"
#include "sca_ref_bias.hpp"
#include <cstdlib>

namespace {
// random symmetric K-bit matrix on the first n spins (same construction as mb::random_k, restricted to i, j < n)
mb::Matrix random_kn(uint64_t seed, double density, int maxabs, int n) {
    mb::Matrix m;
    m.source = "randn:" + std::to_string(seed) + ":" + std::to_string(density) + ":" + std::to_string(maxabs) + ":" + std::to_string(n);
    uint64_t st = seed;
    const uint64_t thr = uint64_t(density * 18446744073709551616.0 - 1.0);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const uint64_t a = sca::Rng::splitmix(st), b = sca::Rng::splitmix(st);
            if (a > thr) continue;
            const int mag = 1 + int(b % uint64_t(maxabs)); const int v = ((b >> 32) & 1) ? mag : -mag;
            mb::at(m, i, j) = int8_t(v); mb::at(m, j, i) = int8_t(v);
            m.edges.push_back({i, j, -v});
        }
    mb::finish(m);
    return m;
}
}  // namespace

int main(int argc, char** argv) {
    using namespace sca;
    if (argc < 4) { std::fprintf(stderr, "usage: gen_vectors_mbb K matrix outdir [suite [case]]\n"); return 2; }
    try {
        const int K = std::atoi(argv[1]);
        const std::string spec = argv[2];
        int n_act = N;
        mb::Matrix mat;
        if (spec.rfind("randn:", 0) == 0) {
            std::stringstream ss(spec.substr(6)); std::string a, b, c, d;
            std::getline(ss, a, ':'); std::getline(ss, b, ':'); std::getline(ss, c, ':'); std::getline(ss, d, ':');
            n_act = std::stoi(d);
            if (n_act < 1 || n_act > N) throw std::runtime_error("randn: n must be 1..N");
            mat = random_kn(std::stoull(a), std::stod(b), std::stoi(c), n_act);
        } else mat = mb::load(spec);
        const std::string dir = argv[3];
        const std::string suite = argc > 4 ? argv[4] : "short";
        const bool k2000 = mat.from_packed;
        const double auto_scale = mb::j2_scale(mat);
        long nnz = 0; for (int8_t v : mat.J) nnz += v != 0;
        const double qd = std::round(2.0 * (0.3 * double(nnz) / N + 0.8)) / 2.0;   // ternary test suites only
        if (K < 2) throw std::runtime_error("K >= 2");
        const std::vector<uint64_t> P = mb::planes(mat, K);
        const Graph g = mb::graph_of(mat);
        FILE* fs = std::fopen((dir + "/stim.hex").c_str(), "w");
        FILE* fe = std::fopen((dir + "/expect.hex").c_str(), "w");
        FILE* fc = std::fopen((dir + "/counts.txt").c_str(), "w");
        FILE* fi = std::fopen((dir + "/cases.txt").c_str(), "w");
        if (!fs || !fe || !fc || !fi) { std::perror("open"); return 1; }
        auto beat = [&](uint64_t hi, uint64_t lo) { std::fprintf(fs, "%016llx%016llx\n", (unsigned long long)hi, (unsigned long long)lo); };
        auto word = [&](uint64_t v) { std::fprintf(fe, "%016llx\n", (unsigned long long)v); };

        // use_b: 0 = no bias (b = 0), 1 = the generated bias vector
        struct Case { const char* name; mb::Sched s; int trials; int trace; int keep; int use_b; };
        auto S = [](double t0, double t1, int steps, double q, double lam, bool ramp, double jv, double kappa, double scale) {
            mb::Sched c; c.t0 = t0; c.t1 = t1; c.S = steps; c.q = q; c.lam = lam; c.ramp = ramp; c.jv = jv; c.kappa = kappa;
            c.scale = scale; return c;
        };
        std::vector<Case> cases;
        const bool bias_suite = suite == "bias";
        if (k2000 && suite == "short")   // exactly the v6.4 suite (gen_vectors.cpp): same tables, seeds and trial ids
            cases = {{"plain_S40", S(30, 5, 40, 4, 0, false, 0, 0, 1), 3, 1, 0, 0},
                     {"onsager_ramp_S40", S(20, 5, 40, 4, 0.7, true, 0, 0, 1), 3, 1, 0, 0},
                     {"tec_S40", S(20, 5, 40, 8, 0, false, -8, 0, 1), 2, 0, 0, 0},
                     {"tec_S40_keep", S(20, 5, 40, 8, 0, false, -8, 0, 1), 2, 1, 1, 0}};
        else if (k2000 && suite == "proto")   // gen_vectors.cpp "proto" (same order, seeds, trial ids), plus X5 as case 6
            cases = {{"O1", S(12, 5, 960, 8, 1.05, true, 0, 0, 1), 1, 0, 0, 0},
                     {"O2", S(15, 5, 560, 6, 0.9, true, 0, 0, 1), 1, 0, 0, 0},
                     {"O3", S(12, 5, 560, 6, 0.9, true, 0, 0, 1), 1, 0, 0, 0},
                     {"O4", S(12, 5, 360, 6, 0.9, false, 0, 0, 1), 1, 0, 0, 0},
                     {"P1", S(30, 5, 1560, 4, 0, false, 0, 0, 1), 1, 0, 0, 0},
                     {"T1", S(30, 5, 1560, 8, 0, false, -4, 0, 1), 1, 0, 0, 0},
                     {"X5", S(15, 5, 280, 8, 0, true, 0, 1.75, 1), 1, 0, 0, 0}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "short")
            cases = {{"plain_S40", S(4, 0.5, 40, qd, 0, false, 0, 0, 1), 3, 1, 0, 0},
                     {"onsager_ramp_S40", S(4, 0.5, 40, qd, 0.9, true, 0, 0, auto_scale), 3, 1, 0, 0},
                     {"tec_S40", S(4, 0.5, 40, qd, 0, false, -1, 0, 1), 2, 0, 0, 0},
                     {"tec_S40_keep", S(4, 0.5, 40, qd, 0, false, -1, 0, 1), 2, 1, 1, 0}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "neg")
            cases = {{"plain_q1_S40", S(3, 0.3, 40, 1, 0, false, 0, 0, 1), 3, 1, 0, 0},
                     {"onsager_q1_S40", S(2, 0.3, 40, 1, 0.9, true, 0, 0, auto_scale), 3, 1, 0, 0},
                     {"tec_q1_S40", S(2, 0.3, 40, 1, 0, false, -0.5, 0, 1), 2, 0, 0, 0}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "g22ons")   // one trial of the V80 MB power protocol's G22_ons
            // configuration (PROTOCOL_HW_MB_power.txt: q 5, lambda 0.5, ramp, t0 5, t1 0.5, 1024 steps, corr-scale auto)
            cases = {{"G22_ons_S1024", S(5, 0.5, 1024, 5, 0.5, true, 0, 0, auto_scale), 1, 0, 0, 0}};
        else if (!k2000 && mat.maxabs <= 1 && suite == "long")
            cases = {{"onsager_ramp_S600", S(4, 0.5, 600, qd, 0.9, true, 0, 0, auto_scale), 2, 1, 0, 0},
                     {"tecT_ramp_S600", S(4, 0.5, 600, qd, 0, true, 0, 1.5, auto_scale), 1, 1, 0, 0}};
        else if (!k2000 && suite == "short")   // multi-bit random matrix (e.g. K = 4, |J| <= 7)
            cases = {{"plain_S40", S(150, 20, 40, 20, 0, false, 0, 0, 1), 3, 1, 0, 0},
                     {"onsager_ramp_S40", S(100, 20, 40, 20, 0.7, true, 0, 0, auto_scale), 3, 1, 0, 0},
                     {"tec_S40", S(100, 20, 40, 20, 0, false, -20, 0, 1), 2, 0, 0, 0},
                     {"tec_S40_keep", S(100, 20, 40, 20, 0, false, -20, 0, 1), 2, 1, 1, 0}};
        else if (!k2000 && suite == "long")
            cases = {{"onsager_ramp_S400", S(150, 10, 400, 10, 0.9, true, 0, 0, auto_scale), 2, 1, 0, 0},
                     {"tecT_ramp_S300", S(150, 10, 300, 20, 0, true, 0, 1.5, auto_scale), 1, 1, 0, 0}};
        else if (!k2000 && bias_suite && mat.maxabs > 1)   // K-bit random matrix with bias (fields ~ sqrt(n d E[J^2]))
            cases = {{"bias_plain_S40", S(150, 20, 40, 20, 0, false, 0, 0, 1), 3, 1, 0, 1},
                     {"bias_onsager_ramp_S60", S(100, 20, 60, 20, 0.7, true, 0, 0, auto_scale), 3, 1, 0, 1},
                     {"nobias_plain_S40", S(150, 20, 40, 20, 0, false, 0, 0, 1), 2, 1, 0, 0},
                     {"bias_tec_S40", S(100, 20, 40, 20, 0, false, -20, 0, 1), 2, 0, 0, 1},
                     {"bias_tec_S40_keep", S(100, 20, 40, 20, 0, false, -20, 0, 1), 2, 1, 1, 1}};
        else if (!k2000 && bias_suite)                     // ternary random matrix with bias
            cases = {{"bias_plain_S40", S(4, 0.5, 40, qd, 0, false, 0, 0, 1), 3, 1, 0, 1},
                     {"bias_onsager_ramp_S60", S(4, 0.5, 60, qd, 0.9, true, 0, 0, auto_scale), 3, 1, 0, 1},
                     {"nobias_plain_S40", S(4, 0.5, 40, qd, 0, false, 0, 0, 1), 2, 1, 0, 0},
                     {"bias_tec_S40", S(4, 0.5, 40, qd, 0, false, -1, 0, 1), 2, 0, 0, 1},
                     {"bias_tec_S40_keep", S(4, 0.5, 40, qd, 0, false, -1, 0, 1), 2, 1, 1, 1}};
        else throw std::runtime_error("unknown suite " + suite);
        if (argc > 5) { Case one = cases.at(std::atoi(argv[5])); one.keep = 0; cases = {one}; }

        // bias vector (bias suite): most |b| comparable to the fields, ~3% large (16384..32767) to exercise the 17-bit
        // fields; padding spins (i >= n) have b = 0. Pure function of the matrix spec (splitmix64 stream).
        std::vector<int32_t> bias(NP, 0);
        int32_t bmax = 0;
        if (bias_suite) {
            uint64_t st = 0xB1A5ull ^ mb::fingerprint(mat);
            const int br = mat.maxabs > 1 ? 300 : 10;
            for (int i = 0; i < n_act; ++i) {
                const uint64_t a = sca::Rng::splitmix(st), b = sca::Rng::splitmix(st);
                int v;
                if (a % 100 < 3) v = 16384 + int(b % 16384u);                       // large: 16384..32767
                else v = int(b % uint64_t(2 * br + 1)) - br;                        // small: -br..br
                if ((a >> 40) & 1) v = -v;
                bias[i] = v; bmax = std::max(bmax, std::abs(v));
            }
            // range check of the device fields (HB = 17 bits): max_y sum_x |J_xy| + max |b| < 2^16
            int64_t rowmax = 0;
            for (int y = 0; y < NP; ++y) { int64_t r = 0; for (int x = 0; x < NP; ++x) r += std::abs(mb::at(mat, x, y)); rowmax = std::max(rowmax, r); }
            if (rowmax + bmax >= 65536) throw std::runtime_error("bias + couplings exceed the 17-bit field range");
        }

        std::fprintf(fi, "K %d matrix %s n %d maxabs %d sumw %lld fingerprint %016llx auto_scale %.12g max|b| %d\n", K,
                     mat.source.c_str(), n_act, mat.maxabs, (long long)mat.sumw, (unsigned long long)mb::fingerprint(mat),
                     auto_scale, bmax);
        long nbeats = 0, nwords = 0;
        for (size_t ci = 0; ci < cases.size(); ++ci) {
            const Case& c = cases[ci];
            const Tables tb = mb::tables(c.s);
            mb::check_tables(tb, K);
            const uint64_t seed = 0xC0FFEEull + ci; const uint32_t off = 1000 + 100 * uint32_t(ci);
            const int loadJ = ci == 0;
            const int loadB = bias_suite && ci == 0;
            const int raw = bias_suite ? 1 : 0;
            const uint64_t nfield = bias_suite ? uint64_t(n_act) : 0;   // 0 = N (the b = 0 suites keep the v6.4 header)
            beat(seed, (uint64_t(raw) << 37) | (uint64_t(c.use_b) << 36) | (uint64_t(loadB) << 35) | (uint64_t(c.keep) << 34) |
                       (uint64_t(c.trace) << 33) | (uint64_t(loadJ) << 32) | (uint64_t(c.s.S) << 16) | uint64_t(c.trials));
            beat(uint64_t(g.sumw), (nfield << 32) | uint64_t(off));
            nbeats += 2;
            if (loadJ) {
                for (size_t r = 0; r < size_t(K) * NP; ++r)
                    for (int cc = 0; cc < 16; ++cc) beat(P[r * WORDS + 2 * cc + 1], P[r * WORDS + 2 * cc]);
                nbeats += long(K) * NP * 16;
            }
            if (loadB) {
                for (int a = 0; a < NP / 8; ++a) {
                    uint64_t lo = 0, hi = 0;
                    for (int i = 0; i < 4; ++i) lo |= uint64_t(uint16_t(int16_t(bias[8 * a + i]))) << (16 * i);
                    for (int i = 0; i < 4; ++i) hi |= uint64_t(uint16_t(int16_t(bias[8 * a + 4 + i]))) << (16 * i);
                    beat(hi, lo);
                }
                nbeats += NP / 8;
            }
            for (int t = 0; t < (c.keep ? 0 : c.s.S); ++t)
                beat((uint64_t(uint32_t(tb.kcorr[t])) << 32) | uint64_t(uint32_t(tb.kconst[t])),
                     (uint64_t(uint32_t(tb.q[t])) << 32) | uint64_t(uint32_t(tb.fourT[t])));
            nbeats += c.keep ? 0 : c.s.S;
            std::fprintf(fi, "case %zu %s S %d t0 %g t1 %g q %g lambda %g ramp %d jv %g kappa %g scale %.12g trials %d trace %d keep %d bias %d raw %d n %d seed %llu offset %u\n",
                         ci, c.name, c.s.S, c.s.t0, c.s.t1, c.s.q, c.s.lam, int(c.s.ramp), c.s.jv, c.s.kappa, c.s.scale, c.trials,
                         c.trace, c.keep, c.use_b, raw, n_act, (unsigned long long)seed, off);
            for (int tr = 0; tr < c.trials; ++tr) {
                if (!bias_suite) {
                    const Result r = run_trial(g, tb, seed, off + tr, true, mat.J.data());
                    for (int w = 0; w < WORDS; ++w) word(r.spins[w]);
                    word(uint64_t(r.cut)); word(uint64_t(r.flips));
                    nwords += WORDS + 2;
                    if (c.trace) { for (int t = 0; t < c.s.S; ++t) word(uint64_t(uint32_t(r.n_lin[t]))); nwords += c.s.S; }
                    const int64_t ind = mb::cut_of(mat, r.spins.data());
                    if (ind != r.cut) throw std::runtime_error("independent cut disagrees with the reference");
                    std::printf("%-22s trial %d: cut %lld flips %lld\n", c.name, tr, (long long)r.cut, (long long)r.flips);
                    std::fprintf(fi, "  trial %u cut %lld flips %lld\n", off + tr, (long long)r.cut, (long long)r.flips);
                } else {
                    const ResultBias r = run_trial_bias(n_act, mat.J.data(), c.use_b ? bias.data() : nullptr, tb, seed, off + tr, true);
                    for (int w = 0; w < WORDS; ++w) word(r.spins[w]);
                    word(uint64_t(r.sum_sh)); word(uint64_t(r.flips));
                    nwords += WORDS + 2;
                    if (c.trace) { for (int t = 0; t < c.s.S; ++t) word(uint64_t(uint32_t(r.n_lin[t]))); nwords += c.s.S; }
                    // independent check: sum_i s_i h_i from scratch with the biased fields
                    int64_t sh = 0;
                    for (int i = 0; i < n_act; ++i) {
                        const int si = ((r.spins[i >> 6] >> (i & 63)) & 1) ? 1 : -1;
                        int64_t hi = c.use_b ? bias[i] : 0;
                        for (int j = 0; j < n_act; ++j) if (j != i) hi += int64_t(mb::at(mat, i, j)) * (((r.spins[j >> 6] >> (j & 63)) & 1) ? 1 : -1);
                        sh += si * hi;
                    }
                    if (sh != r.sum_sh) throw std::runtime_error("independent sum_sh disagrees with the reference");
                    std::printf("%-22s trial %d: sum_sh %lld sum_bs %lld flips %lld\n", c.name, tr, (long long)r.sum_sh,
                                (long long)r.sum_bs, (long long)r.flips);
                    std::fprintf(fi, "  trial %u sum_sh %lld sum_bs %lld flips %lld\n", off + tr, (long long)r.sum_sh,
                                 (long long)r.sum_bs, (long long)r.flips);
                }
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
