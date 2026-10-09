// Ising instance generator and checker for the bias study (writes SCAJINT8 couplings + SCABIAS1 biases + JSON metadata).
// Model (MULTIBIT_SPEC "Bias"): H(s) = -sum_{i<j} J_ij s_i s_j - sum_i b_i s_i, s in {-1,+1}, x = (1 + s)/2.
//   gen_ising random --n N --kbits K --pzero P --bmax B [--pbig Q --bsmall b] --seed S --out PREFIX
//       symmetric J (P(J=0) = P, else |J| uniform in 1..2^(K-1)-1, random sign); b uniform in [-B, B] with probability Q
//       (default 1), otherwise uniform in [-b, b] (splitmix64 streams)
//   gen_ising zero --n N --out FILE.bias                     all-zero bias file
//   gen_ising gp --n N --deg D --k K --A a --B b --seed S --out PREFIX [--brute 1]
//       graph partitioning (Lucas 2014, unbalanced sizes): random D-regular graph, minimise the cut subject to |{x = 1}| = K,
//       f(x) = b * cut(x) + a * (sum x - K)^2;  J_ij = -a + b [ij in E], b_i = a (2K - n);  H = 2 f + C
//   gen_ising tsp --n N --wmax W --A a --B b --seed S --out PREFIX [--brute 1]
//       TSP (Lucas 2014), city 0 fixed at position 0, x_{v,p} for v, p in 1..n-1 (spin (v-1)(n-1) + (p-1)), symmetric
//       distances uniform in 1..W; f = a sum_v (1 - sum_p x_vp)^2 + a sum_p (1 - sum_v x_vp)^2 + b * tour length;  H = 4 f + C
//   gen_ising brute --jint8 F [--bias F]                    exhaustive ground state (n <= 30, Gray code)
// --brute 1 additionally checks H = sigma f + C on random states, enumerates all 2^n spin states for the Ising minimum and
// checks it against a direct enumeration of the original problem (all K-subsets / all tours).
#include "jmat.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

static uint64_t splitmix(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

struct Inst {
    jmat::IntJ J;
    std::vector<int32_t> b = std::vector<int32_t>(sca::NP, 0);
};

static int64_t ising_energy(const Inst& in, const std::vector<int>& s) {
    int64_t e = 0;
    for (int i = 0; i < in.J.N; ++i) {
        e -= int64_t(in.b[i]) * s[i];
        for (int j = i + 1; j < in.J.N; ++j) e -= int64_t(in.J.at(i, j)) * s[i] * s[j];
    }
    return e;
}

// exhaustive minimum over all 2^n states (Gray code with incremental fields); returns min energy and up to 64 ground states
struct Brute { int64_t emin = INT64_MAX; int64_t count = 0; std::vector<std::vector<int>> gs; };
static Brute brute_force(const Inst& in) {
    const int n = in.J.N;
    if (n > 30) throw std::runtime_error("brute force limited to n <= 30");
    std::vector<int> s(n, -1);
    std::vector<int64_t> h(n);
    for (int i = 0; i < n; ++i) { h[i] = in.b[i]; for (int j = 0; j < n; ++j) if (j != i) h[i] += int64_t(in.J.at(i, j)) * s[j]; }
    int64_t e = ising_energy(in, s);
    Brute br;
    auto visit = [&]() {
        if (e < br.emin) { br.emin = e; br.count = 0; br.gs.clear(); }
        if (e == br.emin) { ++br.count; if (br.gs.size() < 64) br.gs.push_back(s); }
    };
    visit();
    for (uint64_t g = 1; g < (uint64_t(1) << n); ++g) {
        const int k = __builtin_ctzll(g);           // flip spin k: dE = 2 s_k h_k
        e += 2 * int64_t(s[k]) * h[k];
        s[k] = -s[k];
        for (int j = 0; j < n; ++j) if (j != k) h[j] += 2 * int64_t(in.J.at(j, k)) * s[k];
        visit();
    }
    return br;
}

static void write_inst(const Inst& in, const std::string& prefix, const std::string& meta) {
    jmat::write_jint8(in.J, prefix + ".jint8");
    jmat::write_bias(in.b, in.J.N, prefix + ".bias");
    FILE* f = std::fopen((prefix + ".json").c_str(), "w");
    if (!f) throw std::runtime_error("cannot create " + prefix + ".json");
    std::fputs(meta.c_str(), f); std::fputs("\n", f); std::fclose(f);
    std::cout << meta << "\n";
}

static void set_pair(Inst& in, int i, int j, int v) {
    if (v < -127 || v > 127) throw std::runtime_error("coupling outside int8");
    in.J.dense[size_t(i) * sca::NP + j] = in.J.dense[size_t(j) * sca::NP + i] = int8_t(v);
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) throw std::runtime_error("usage: gen_ising random|zero|gp|tsp|brute ...");
        const std::string cmd = argv[1];
        std::string out, jpath, bpath; int n = 0, K = 2, deg = 3, kk = 0, wmax = 2, brute = 0; double pzero = 0; int64_t bmax = 0, A = 1, Bw = 1;
        uint64_t seed = 1; double pbig = 1.0; int64_t bsmall = 0;
        for (int a = 2; a < argc; ++a) {
            const std::string k = argv[a];
            auto val = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++a]); };
            if (k == "--out") out = val(); else if (k == "--n") n = std::stoi(val()); else if (k == "--kbits") K = std::stoi(val());
            else if (k == "--pzero") pzero = std::stod(val()); else if (k == "--bmax") bmax = std::stoll(val());
            else if (k == "--seed") seed = std::stoull(val()); else if (k == "--deg") deg = std::stoi(val());
            else if (k == "--k") kk = std::stoi(val()); else if (k == "--A") A = std::stoll(val()); else if (k == "--B") Bw = std::stoll(val());
            else if (k == "--wmax") wmax = std::stoi(val()); else if (k == "--brute") brute = std::stoi(val());
            else if (k == "--jint8") jpath = val(); else if (k == "--bias") bpath = val();
            else if (k == "--pbig") pbig = std::stod(val()); else if (k == "--bsmall") bsmall = std::stoll(val());
            else throw std::runtime_error("unknown argument " + k);
        }
        std::ostringstream meta;
        if (cmd == "zero") {
            if (out.empty() || n < 1 || n > sca::NP) throw std::runtime_error("zero: --n, --out");
            jmat::write_bias(std::vector<int32_t>(sca::NP, 0), n, out);
            std::cout << "{\"wrote\":\"" << out << "\",\"N\":" << n << "}\n";
            return 0;
        }
        if (cmd == "brute") {
            Inst in; in.J = jmat::from_jint8(jpath);
            if (!bpath.empty()) in.b = jmat::load_bias(bpath, in.J.N);
            const Brute br = brute_force(in);
            std::cout << "{\"N\":" << in.J.N << ",\"emin\":" << br.emin << ",\"ground_states\":" << br.count << "}\n";
            return 0;
        }
        if (cmd == "random") {
            if (out.empty() || K < 2 || K > 8 || n < 1 || n > sca::NP) throw std::runtime_error("random: --out, 2<=K<=8, 1<=n<=2048");
            const int M = (1 << (K - 1)) - 1;
            Inst in; in.J.N = n; in.J.source = "random";
            uint64_t st = seed ^ 0xA5A5F00D12345678ull, sb = seed ^ 0x5B1A5B1A5B1A5B1Aull;
            const double lim = std::ldexp(pzero, 64);
            const uint64_t zthr = lim >= 18446744073709551615.0 ? ~0ull : uint64_t(lim);
            for (int x = 0; x < n; ++x)
                for (int y = x + 1; y < n; ++y) {
                    const uint64_t r0 = splitmix(st), r1 = splitmix(st);
                    int v = 0;
                    if (!(pzero > 0 && r0 < zthr)) { const int mag = 1 + int((r1 >> 1) % uint64_t(M)); v = (r1 & 1) ? -mag : mag; }
                    set_pair(in, x, y, v);
                }
            for (int i = 0; i < n; ++i) {
                const bool big = pbig >= 1.0 || double(splitmix(sb) >> 11) * 0x1.0p-53 < pbig;
                const int64_t bm = big ? bmax : bsmall;
                in.b[i] = bm > 0 ? int32_t(int64_t(splitmix(sb) % uint64_t(2 * bm + 1)) - bm) : 0;
            }
            jmat::finish(in.J);
            int64_t bm = 0; for (int i = 0; i < n; ++i) bm = std::max<int64_t>(bm, std::abs(in.b[i]));
            meta << "{\"type\":\"random\",\"N\":" << n << ",\"kbits\":" << K << ",\"pzero\":" << pzero << ",\"bmax_param\":" << bmax << ",\"pbig\":" << pbig
                 << ",\"bsmall\":" << bsmall << ",\"seed\":" << seed
                 << ",\"maxabs\":" << in.J.maxabs << ",\"hmax\":" << in.J.hmax << ",\"bmax\":" << bm << ",\"sumw\":" << in.J.sumw << "}";
            write_inst(in, out, meta.str());
            return 0;
        }
        if (cmd == "gp") {
            if (out.empty() || n < 4 || n > sca::NP || kk < 1 || kk >= n || (n * deg) % 2) throw std::runtime_error("gp: bad arguments");
            // random deg-regular simple graph by the pairing model with rejection
            std::set<std::pair<int, int>> E; uint64_t st = seed ^ 0x6A09E667F3BCC908ull; int attempts = 0;
            for (;;) {
                ++attempts; E.clear();
                std::vector<int> stubs; for (int v = 0; v < n; ++v) for (int d = 0; d < deg; ++d) stubs.push_back(v);
                for (int i = int(stubs.size()) - 1; i > 0; --i) std::swap(stubs[i], stubs[splitmix(st) % uint64_t(i + 1)]);
                bool ok = true;
                for (size_t i = 0; i < stubs.size(); i += 2) {
                    int u = stubs[i], v = stubs[i + 1]; if (u > v) std::swap(u, v);
                    if (u == v || E.count({u, v})) { ok = false; break; }
                    E.insert({u, v});
                }
                if (ok) break;
                if (attempts > 100000) throw std::runtime_error("no simple regular graph");
            }
            Inst in; in.J.N = n; in.J.source = "gp";
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) set_pair(in, i, j, int(-A + Bw * (E.count({i, j}) ? 1 : 0)));
            for (int i = 0; i < n; ++i) in.b[i] = int32_t(A * (2 * kk - n));
            jmat::finish(in.J);
            auto f_of = [&](const std::vector<int>& s) {   // f(x) = B cut + A (sum x - k)^2
                int64_t cut = 0, sx = 0;
                for (auto& e : E) cut += s[e.first] != s[e.second];
                for (int i = 0; i < n; ++i) sx += s[i] > 0;
                return std::make_pair(Bw * cut + A * (sx - kk) * (sx - kk), std::make_pair(cut, sx));
            };
            std::vector<int> s0(n, -1);
            const int64_t C = ising_energy(in, s0) - 2 * f_of(s0).first;
            uint64_t sr = seed ^ 0x1234; int bad = 0;
            for (int t = 0; t < 200; ++t) { std::vector<int> s(n); for (auto& v : s) v = (splitmix(sr) & 1) ? 1 : -1; bad += ising_energy(in, s) != 2 * f_of(s).first + C; }
            if (bad) throw std::runtime_error("gp: H != 2f + C on random states");
            meta << "{\"type\":\"graph_partitioning\",\"N\":" << n << ",\"deg\":" << deg << ",\"k\":" << kk << ",\"A\":" << A << ",\"B\":" << Bw << ",\"seed\":" << seed
                 << ",\"edges\":" << E.size() << ",\"sigma\":2,\"C\":" << C << ",\"relation\":\"H = 2 f + C, f = B*cut + A*(sum x - k)^2, x_i = (1+s_i)/2\""
                 << ",\"maxabs\":" << in.J.maxabs << ",\"hmax\":" << in.J.hmax << ",\"bmax\":" << std::llabs(A * (2 * kk - n));
            if (brute) {
                if (n > 30) throw std::runtime_error("gp brute force limited to n <= 30");
                const Brute br = brute_force(in);
                int64_t best = INT64_MAX, nbest = 0;   // direct: all k-subsets
                for (uint64_t m = 0; m < (uint64_t(1) << n); ++m) {
                    if (__builtin_popcountll(m) != kk) continue;
                    int64_t cut = 0; for (auto& e : E) cut += ((m >> e.first) & 1) != ((m >> e.second) & 1);
                    if (cut < best) { best = cut; nbest = 0; }
                    nbest += cut == best;
                }
                bool all_feasible_opt = true;
                for (auto& g : br.gs) { auto fv = f_of(g); all_feasible_opt &= fv.second.second == kk && fv.second.first == best; }
                const bool ok = br.emin == 2 * Bw * best + C && all_feasible_opt && br.count == nbest;
                meta << ",\"brute\":{\"ising_emin\":" << br.emin << ",\"ising_ground_states\":" << br.count << ",\"direct_min_cut\":" << best
                     << ",\"direct_optima\":" << nbest << ",\"emin_expected\":" << 2 * Bw * best + C << ",\"ground_states_feasible_optimal\":" << (all_feasible_opt ? "true" : "false")
                     << ",\"consistent\":" << (ok ? "true" : "false") << "}";
                if (!ok) { meta << "}"; std::cout << meta.str() << "\n"; throw std::runtime_error("gp brute-force check failed"); }
            }
            meta << "}";
            write_inst(in, out, meta.str());
            FILE* fe = std::fopen((out + ".edges").c_str(), "w");
            for (auto& e : E) std::fprintf(fe, "%d %d\n", e.first, e.second);
            std::fclose(fe);
            return 0;
        }
        if (cmd == "tsp") {
            if (out.empty() || n < 3 || (n - 1) * (n - 1) > sca::NP || wmax < 1) throw std::runtime_error("tsp: bad arguments");
            const int m = n - 1, nv = m * m;
            uint64_t st = seed ^ 0x3C6EF372FE94F82Bull;
            std::vector<std::vector<int>> W(n, std::vector<int>(n, 0));
            for (int u = 0; u < n; ++u) for (int v = u + 1; v < n; ++v) W[u][v] = W[v][u] = 1 + int(splitmix(st) % uint64_t(wmax));
            auto idx = [&](int v, int p) { return (v - 1) * m + (p - 1); };
            std::vector<std::vector<int64_t>> Q(nv, std::vector<int64_t>(nv, 0)); std::vector<int64_t> c(nv, 0); int64_t c0 = 0;
            auto addQ = [&](int i, int j, int64_t v) { if (i > j) std::swap(i, j); Q[i][j] += v; };
            for (int v = 1; v <= m; ++v) {   // row constraint (city v visited once) and column constraint (position v filled once)
                c0 += 2 * A;
                for (int p = 1; p <= m; ++p) {
                    c[idx(v, p)] -= A; c[idx(p, v)] -= A;
                    for (int p2 = p + 1; p2 <= m; ++p2) { addQ(idx(v, p), idx(v, p2), 2 * A); addQ(idx(p, v), idx(p2, v), 2 * A); }
                }
            }
            for (int p = 1; p + 1 <= m; ++p)
                for (int u = 1; u <= m; ++u) for (int v = 1; v <= m; ++v) if (u != v) addQ(idx(u, p), idx(v, p + 1), Bw * W[u][v]);
            for (int v = 1; v <= m; ++v) { c[idx(v, 1)] += Bw * W[0][v]; c[idx(v, m)] += Bw * W[v][0]; }
            // Ising with sigma = 4: H = 4 f + C, J_ij = -Q_ij, b_i = -(2 c_i + sum_j Q_ij)
            Inst in; in.J.N = nv; in.J.source = "tsp";
            for (int i = 0; i < nv; ++i) {
                int64_t rs = 0;
                for (int j = 0; j < nv; ++j) if (j != i) rs += i < j ? Q[i][j] : Q[j][i];
                for (int j = i + 1; j < nv; ++j) set_pair(in, i, j, int(-Q[i][j]));
                const int64_t bi = -(2 * c[i] + rs);
                if (std::llabs(bi) > 32767) throw std::runtime_error("tsp bias outside 16 bits");
                in.b[i] = int32_t(bi);
            }
            jmat::finish(in.J);
            auto f_of = [&](const std::vector<int>& s) {
                int64_t f = c0;
                for (int i = 0; i < nv; ++i) if (s[i] > 0) { f += c[i]; for (int j = i + 1; j < nv; ++j) if (s[j] > 0) f += Q[i][j]; }
                return f;
            };
            std::vector<int> s0(nv, -1);
            const int64_t C = ising_energy(in, s0) - 4 * f_of(s0);
            uint64_t sr = seed ^ 0x777; int bad = 0;
            for (int t = 0; t < 200; ++t) { std::vector<int> s(nv); for (auto& v : s) v = (splitmix(sr) & 1) ? 1 : -1; bad += ising_energy(in, s) != 4 * f_of(s) + C; }
            if (bad) throw std::runtime_error("tsp: H != 4f + C on random states");
            // decode helper: valid permutation -> tour length, else -1
            auto tour_len = [&](const std::vector<int>& s) -> int64_t {
                std::vector<int> at(m + 1, -1);
                for (int v = 1; v <= m; ++v) for (int p = 1; p <= m; ++p) if (s[idx(v, p)] > 0) { if (at[p] >= 0) return -1; at[p] = v; }
                std::vector<int> seen(m + 1, 0);
                for (int p = 1; p <= m; ++p) { if (at[p] < 0 || seen[at[p]]) return -1; seen[at[p]] = 1; }
                int64_t L = W[0][at[1]] + W[at[m]][0];
                for (int p = 1; p < m; ++p) L += W[at[p]][at[p + 1]];
                return L;
            };
            meta << "{\"type\":\"tsp\",\"cities\":" << n << ",\"N\":" << nv << ",\"wmax\":" << wmax << ",\"A\":" << A << ",\"B\":" << Bw << ",\"seed\":" << seed
                 << ",\"sigma\":4,\"C\":" << C << ",\"relation\":\"H = 4 f + C; valid tour: H = 4*B*length + C; spin (v-1)(n-1)+(p-1) = city v at position p; city 0 at position 0\""
                 << ",\"maxabs\":" << in.J.maxabs << ",\"hmax\":" << in.J.hmax;
            int64_t bm = 0; for (int i = 0; i < nv; ++i) bm = std::max<int64_t>(bm, std::abs(in.b[i]));
            meta << ",\"bmax\":" << bm << ",\"W\":[";
            for (int u = 0; u < n; ++u) { meta << (u ? "," : "") << "["; for (int v = 0; v < n; ++v) meta << (v ? "," : "") << W[u][v]; meta << "]"; }
            meta << "]";
            if (brute) {
                if (nv > 30) throw std::runtime_error("tsp brute force limited to (n-1)^2 <= 30");
                const Brute br = brute_force(in);
                std::vector<int> perm(m); std::iota(perm.begin(), perm.end(), 1);
                int64_t best = INT64_MAX, nbest = 0;
                do {
                    int64_t L = W[0][perm[0]] + W[perm[m - 1]][0];
                    for (int p = 0; p + 1 < m; ++p) L += W[perm[p]][perm[p + 1]];
                    if (L < best) { best = L; nbest = 0; }
                    nbest += L == best;
                } while (std::next_permutation(perm.begin(), perm.end()));
                bool all_opt = true; for (auto& g : br.gs) all_opt &= tour_len(g) == best;
                const bool ok = br.emin == 4 * Bw * best + C && all_opt && br.count == nbest;
                meta << ",\"brute\":{\"ising_emin\":" << br.emin << ",\"ising_ground_states\":" << br.count << ",\"direct_min_length\":" << best
                     << ",\"direct_optimal_tours\":" << nbest << ",\"emin_expected\":" << 4 * Bw * best + C << ",\"ground_states_valid_optimal\":" << (all_opt ? "true" : "false")
                     << ",\"consistent\":" << (ok ? "true" : "false") << "}";
                if (!ok) { meta << "}"; std::cout << meta.str() << "\n"; throw std::runtime_error("tsp brute-force check failed"); }
            }
            meta << "}";
            write_inst(in, out, meta.str());
            return 0;
        }
        throw std::runtime_error("unknown command " + cmd);
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << "\n"; return 1; }
}
