// Small graph-partitioning and TSP Ising instances with integer couplings and biases (MULTIBIT_SPEC.md, "Bias"), checked
// by brute force: the Ising ground state over all 2^n spin states must be a feasible solution with the problem optimum
// (found by direct enumeration of partitions / tours). H(s) = -sum_{i<j} J_ij s_i s_j - sum_i b_i s_i (+ constant).
// QUBO -> Ising, scaled by 4: x = (s + 1) / 2, H_qubo = sum_{i<j} Q_ij x_i x_j + sum_i c_i x_i gives J_ij = -Q_ij and
// b_i = -(sum_j Q_ij + 2 c_i).
// Build: c++ -O3 -std=c++17 gen_gp_tsp.cpp -o gen_gp_tsp
// Usage: gen_gp_tsp gp <n> <p> <m> <seed> <A> <B> <out>    graph partition, sizes (n + m)/2 and (n - m)/2:
//            H = A (sum s - m)^2 + B sum_E (1 - s_i s_j) / 2  ->  J_ij = (B/2)[ij in E] - 2A, b_i = 2 A m   (B even)
//        gen_gp_tsp tsp <c> <dmax> <seed> <A> <B> <out>    TSP on c cities, x_{v,t} -> spin v c + t:
//            H = A sum_v (1 - sum_t x_vt)^2 + A sum_t (1 - sum_v x_vt)^2 + B sum_{u != v} d_uv sum_t x_ut x_{v,t+1}
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <string>
#include <vector>

struct Ising { int n; std::vector<int> J; std::vector<long> b; };   // J[n*n]

static long energy2(const Ising& is, uint64_t st) {   // 2 H(s) = -2 sum_{i<j} J s s - 2 sum b s  (integer)
    long e = 0;
    for (int i = 0; i < is.n; ++i) {
        const int si = (st >> i) & 1 ? 1 : -1;
        e -= 2 * is.b[i] * si;
        for (int j = i + 1; j < is.n; ++j) { const int sj = (st >> j) & 1 ? 1 : -1; e -= 2L * is.J[i * is.n + j] * si * sj; }
    }
    return e;
}
static uint64_t ground(const Ising& is, long& emin, int& nmin) {
    emin = 1L << 62; nmin = 0; uint64_t arg = 0;
    for (uint64_t st = 0; st < (uint64_t(1) << is.n); ++st) {
        const long e = energy2(is, st);
        if (e < emin) { emin = e; arg = st; nmin = 1; } else if (e == emin) ++nmin;
    }
    return arg;
}
static void write(const Ising& is, const std::string& out, const std::string& hdr, long emin) {
    FILE* f = std::fopen(out.c_str(), "w");
    std::fprintf(f, "# %s\n# brute-force ground state: 2H = %ld (H = %ld)\nn %d\n", hdr.c_str(), emin, emin / 2, is.n);
    for (int i = 0; i < is.n; ++i) for (int j = i + 1; j < is.n; ++j) if (is.J[i * is.n + j]) std::fprintf(f, "J %d %d %d\n", i, j, is.J[i * is.n + j]);
    for (int i = 0; i < is.n; ++i) if (is.b[i]) std::fprintf(f, "b %d %ld\n", i, is.b[i]);
    std::fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    const std::string kind = argv[1];
    if (kind == "gp") {
        const int n = std::atoi(argv[2]); const double p = std::atof(argv[3]); const int m = std::atoi(argv[4]);
        std::mt19937_64 rg(std::strtoull(argv[5], 0, 0)); const int A = std::atoi(argv[6]), B = std::atoi(argv[7]);
        std::uniform_real_distribution<double> U(0, 1);
        std::vector<int> E(n * n, 0);
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (U(rg) < p) E[i * n + j] = E[j * n + i] = 1;
        Ising is{n, std::vector<int>(n * n, 0), std::vector<long>(n, 0)};
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (i != j) is.J[i * n + j] = (B / 2) * E[i * n + j] - 2 * A;
        for (int i = 0; i < n; ++i) is.b[i] = 2L * A * m;
        long emin; int nmin; const uint64_t g = ground(is, emin, nmin);
        // direct optimum: among partitions with sum s = m (part of +1 spins has (n + m) / 2 vertices), minimum cut
        long best = 1L << 40; int nbest = 0;
        for (uint64_t st = 0; st < (uint64_t(1) << n); ++st) {
            if (2 * __builtin_popcountll(st) - n != m) continue;
            long cut = 0;
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (E[i * n + j] && (((st >> i) ^ (st >> j)) & 1)) ++cut;
            if (cut < best) { best = cut; nbest = 1; } else if (cut == best) ++nbest;
        }
        long gcut = 0; const int gsum = 2 * __builtin_popcountll(g) - n;
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (E[i * n + j] && (((g >> i) ^ (g >> j)) & 1)) ++gcut;
        long edges = 0; for (int i = 0; i < n * n; ++i) edges += E[i]; edges /= 2;
        const bool ok = gsum == m && gcut == best;
        char hdr[512];
        std::snprintf(hdr, sizeof hdr, "graph partition n=%d p=%g edges=%ld target sum s=%d A=%d B=%d; optimum cut %ld (%d optimal partitions); "
                      "Ising ground state: sum s=%d cut %ld (%d degenerate) -> %s", n, p, edges, m, A, B, best, nbest, gsum, gcut, nmin,
                      ok ? "MATCHES the partition optimum" : "DOES NOT match");
        write(is, argv[8], hdr, emin);
        std::printf("%s\n", hdr);
        return ok ? 0 : 3;
    }
    if (kind == "tsp") {
        const int c = std::atoi(argv[2]), dmax = std::atoi(argv[3]); std::mt19937_64 rg(std::strtoull(argv[4], 0, 0));
        const int A = std::atoi(argv[5]), B = std::atoi(argv[6]); const int n = c * c;
        std::uniform_int_distribution<int> ud(1, dmax);
        std::vector<int> d(c * c, 0);
        for (int u = 0; u < c; ++u) for (int v = u + 1; v < c; ++v) d[u * c + v] = d[v * c + u] = ud(rg);
        std::vector<long> Q(n * n, 0), cc(n, -2L * A);
        auto id = [&](int v, int t) { return v * c + t; };
        for (int v = 0; v < c; ++v) for (int t = 0; t < c; ++t) for (int t2 = t + 1; t2 < c; ++t2) { Q[id(v, t) * n + id(v, t2)] += 2 * A; Q[id(v, t2) * n + id(v, t)] += 2 * A; }
        for (int t = 0; t < c; ++t) for (int u = 0; u < c; ++u) for (int v = u + 1; v < c; ++v) { Q[id(u, t) * n + id(v, t)] += 2 * A; Q[id(v, t) * n + id(u, t)] += 2 * A; }
        for (int t = 0; t < c; ++t) for (int u = 0; u < c; ++u) for (int v = 0; v < c; ++v) if (u != v) {
            const int a = id(u, t), b = id(v, (t + 1) % c); Q[a * n + b] += B * d[u * c + v]; Q[b * n + a] += B * d[u * c + v]; }
        Ising is{n, std::vector<int>(n * n, 0), std::vector<long>(n, 0)};
        for (int i = 0; i < n; ++i) { long rs = 0; for (int j = 0; j < n; ++j) if (j != i) { is.J[i * n + j] = int(-Q[i * n + j]); rs += Q[i * n + j]; } is.b[i] = -(rs + 2 * cc[i]); }
        // direct optimum over tours (fix city 0 at time 0)
        std::vector<int> perm(c); std::iota(perm.begin(), perm.end(), 0); long best = 1L << 40; int nbest = 0;
        do { long L = 0; for (int t = 0; t < c; ++t) L += d[perm[t] * c + perm[(t + 1) % c]]; if (L < best) { best = L; nbest = 1; } else if (L == best) ++nbest; }
        while (std::next_permutation(perm.begin() + 1, perm.end()));
        long emin; int nmin; uint64_t g = 0;
        if (n <= 30) g = ground(is, emin, nmin); else { emin = 0; nmin = 0; }
        // decode the ground state as a tour
        bool valid = true; std::vector<int> at(c, -1);
        for (int t = 0; t < c; ++t) { int cnt = 0; for (int v = 0; v < c; ++v) if ((g >> id(v, t)) & 1) { ++cnt; at[t] = v; } if (cnt != 1) valid = false; }
        for (int v = 0; v < c && valid; ++v) { int cnt = 0; for (int t = 0; t < c; ++t) if ((g >> id(v, t)) & 1) ++cnt; if (cnt != 1) valid = false; }
        long gl = -1; if (valid) { gl = 0; for (int t = 0; t < c; ++t) gl += d[at[t] * c + at[(t + 1) % c]]; }
        const bool ok = valid && gl == best;
        int jmax = 0; long bmax = 0; for (int i = 0; i < n * n; ++i) jmax = std::max(jmax, std::abs(is.J[i])); for (long v : is.b) bmax = std::max(bmax, std::labs(v));
        char hdr[512];
        std::snprintf(hdr, sizeof hdr, "TSP c=%d cities (n=%d spins) dmax=%d A=%d B=%d max|J|=%d max|b|=%ld; optimal tour length %ld (%d tours, city 0 fixed); "
                      "Ising ground state: %s tour length %ld (%d degenerate states) -> %s", c, n, dmax, A, B, jmax, bmax, best, nbest,
                      valid ? "valid" : "INVALID", gl, nmin, ok ? "MATCHES the tour optimum" : "DOES NOT match");
        write(is, argv[7], hdr, emin);
        std::printf("%s\n", hdr);
        return ok ? 0 : 3;
    }
    return 2;
}
