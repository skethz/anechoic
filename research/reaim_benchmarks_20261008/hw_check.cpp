// Bit-exactness driver: runs the golden reference sca::run_trial_bias (fpga/v80_sca/src/sca_ref_bias.hpp, unedited) on an
// exported instance (SCAJINT8 couplings via gpu/multibit_20261007/src/jmat.hpp, int32 bias) with tables built by
// fpga/v80_sca/src/mb_common.hpp::tables from a schedule, and prints one JSON line per trial:
// flips, sum_sh, sum_bs, the final spins (hex, bit i = 1 means +1) and the n_lin trace (comma list).
// Build (gpu-host): c++ -O2 -std=c++17 -DSCA_LANES=256 -I<fpga/v80_sca/src> -I<gpu/multibit_20261007/src> hw_check.cpp -o hw_check
// Usage: hw_check J.jint8 bias.bin t0 t1 S q lam jv kappa ramp(0/1) seed trial0 ntrials
#include "sca_ref_bias.hpp"
#include "mb_common.hpp"
#include "jmat.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 14) { std::fprintf(stderr, "usage: %s J.jint8 bias.bin t0 t1 S q lam jv kappa ramp seed trial0 ntrials\n", argv[0]); return 2; }
    const jmat::IntJ m = jmat::from_jint8(argv[1]);
    std::vector<int32_t> bias(sca::NP, 0);
    {
        FILE* f = std::fopen(argv[2], "rb");
        if (!f) { std::fprintf(stderr, "cannot open bias\n"); return 2; }
        const size_t got = std::fread(bias.data(), 4, size_t(m.N), f);
        std::fclose(f);
        if (got != size_t(m.N)) { std::fprintf(stderr, "short bias file\n"); return 2; }
    }
    mb::Sched sc;
    sc.t0 = std::atof(argv[3]); sc.t1 = std::atof(argv[4]); sc.S = std::atoi(argv[5]); sc.q = std::atof(argv[6]);
    sc.lam = std::atof(argv[7]); sc.jv = std::atof(argv[8]); sc.kappa = std::atof(argv[9]); sc.ramp = std::atoi(argv[10]) != 0;
    const uint64_t seed = std::strtoull(argv[11], nullptr, 10);
    const uint32_t trial0 = uint32_t(std::atoi(argv[12])), nt = uint32_t(std::atoi(argv[13]));
    const sca::Tables tb = mb::tables(sc);
    // tables (for the Python comparison)
    std::printf("{\"tables\":{\"fourT\":[");
    for (size_t t = 0; t < tb.fourT.size(); ++t) std::printf("%s%d", t ? "," : "", tb.fourT[t]);
    std::printf("],\"q\":[");
    for (size_t t = 0; t < tb.q.size(); ++t) std::printf("%s%d", t ? "," : "", tb.q[t]);
    std::printf("],\"kcorr\":[");
    for (size_t t = 0; t < tb.kcorr.size(); ++t) std::printf("%s%lld", t ? "," : "", (long long)tb.kcorr[t]);
    std::printf("],\"kconst\":[");
    for (size_t t = 0; t < tb.kconst.size(); ++t) std::printf("%s%d", t ? "," : "", tb.kconst[t]);
    std::printf("]}}\n");
    for (uint32_t k = 0; k < nt; ++k) {
        const uint32_t trial = trial0 + k;
        const sca::ResultBias r = sca::run_trial_bias(m.N, m.dense.data(), bias.data(), tb, seed, trial, true);
        std::printf("{\"trial\":%u,\"flips\":%lld,\"sum_sh\":%lld,\"sum_bs\":%lld,\"spins\":\"", trial, (long long)r.flips,
                    (long long)r.sum_sh, (long long)r.sum_bs);
        for (int w = 0; w < sca::WORDS; ++w) std::printf("%016llx", (unsigned long long)r.spins[w]);
        std::printf("\",\"nlin\":[");
        for (size_t t = 0; t < r.n_lin.size(); ++t) std::printf("%s%d", t ? "," : "", r.n_lin[t]);
        std::printf("]}\n");
    }
    return 0;
}
