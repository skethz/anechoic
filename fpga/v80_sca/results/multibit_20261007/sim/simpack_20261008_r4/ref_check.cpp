// Quick check: run_trial_bias on an exported jint8 problem with a given schedule, seed and trials; print sum s, energy,
// flips and, for a square n (TSP one-hot), whether the state is a permutation matrix.
#include "mb_common4.hpp"
#include <cmath>
#include <cstdlib>
int main(int argc, char** argv) {
    // ref_check <jint8 spec> t0 t1 S q lam jv kappa ramp seed trial0 n
    const mb4::Matrix m = mb4::load(argv[1]);
    mb4::Sched sc; sc.t0 = atof(argv[2]); sc.t1 = atof(argv[3]); sc.S = atoi(argv[4]); sc.q = atof(argv[5]); sc.lam = atof(argv[6]);
    sc.jv = atof(argv[7]); sc.kappa = atof(argv[8]); sc.ramp = atoi(argv[9]) != 0;
    const sca::Tables tb = mb4::tables(sc);
    const uint64_t seed = strtoull(argv[10], 0, 10); const int t0 = atoi(argv[11]), nt = atoi(argv[12]);
    const int c = int(std::lround(std::sqrt(double(m.n))));
    int feas = 0;
    for (int k = 0; k < nt; ++k) {
        const sca::ResultBias r = sca::run_trial_bias(m.n, m.J.data(), m.bias.data(), tb, seed, uint32_t(t0 + k), false);
        auto bit = [&](int i) { return int((r.spins[i >> 6] >> (i & 63)) & 1); };
        int ss = 0; for (int i = 0; i < m.n; ++i) ss += bit(i) ? 1 : -1;
        bool perm = c * c == m.n;
        if (perm) for (int v = 0; v < c; ++v) { int a = 0, b = 0; for (int j = 0; j < c; ++j) { a += bit(v * c + j); b += bit(j * c + v); } perm &= a == 1 && b == 1; }
        feas += perm;
        printf("trial %d: sum_s %d energy %lld flips %lld%s\n", t0 + k, ss, (long long)(-(r.sum_sh + r.sum_bs) / 2), (long long)r.flips,
               c * c == m.n ? (perm ? " permutation" : " not-permutation") : "");
    }
    printf("permutations %d of %d\n", feas, nt);
}
