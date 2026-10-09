// Bit-exact testbench: HLS kernel vs sca_ref.hpp on K2000 (plain, Onsager+ramp, TEC; short and full schedules).
#include "sca_ref.hpp"
#include <cmath>
#include <cstdlib>

extern "C" void sca_v80(const uint64_t*, const int32_t*, const int32_t*, const int32_t*, const int64_t*, uint64_t*, int32_t*,
                        int, int, uint64_t, uint32_t, int64_t, int, int);

static std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S); for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1)); return T;
}

int main(int argc, char** argv) {
    using namespace sca;
    const char* path = argc > 1 ? argv[1] : "../../v80_snowball/data/K2000.bin";
    const bool full = std::getenv("SCA_FULL") != nullptr;  // also run a full-length schedule (slower in RTL cosim)
    Graph g = load_graph(path);
    std::vector<uint64_t> J64(size_t(NP) * WORDS, 0);
    for (int x = 0; x < N; ++x) for (int w = 0; w < WORDS; ++w) J64[size_t(x) * WORDS + w] = g.J[x][w];
    struct Case { const char* name; double t0; int S; double q; double lam; bool ramp; double jv; int trials; };
    std::vector<Case> cases = {{"plain_S40", 30, 40, 4, 0, false, 0, 3}, {"onsager_ramp_S40", 20, 40, 4, 0.7, true, 0, 3},
                               {"tec_S40", 20, 40, 8, 0, false, -8, 2}};
    if (full) cases.push_back({"onsager_ramp_S560", 30, 560, 4, 0.7, true, 0, 2});
    if (std::getenv("SCA_COSIM")) cases = {{"cosim_plain_S40", 30, 40, 4, 0, false, 0, 1}, {"cosim_onsager_ramp_S40", 20, 40, 4, 0.7, true, 0, 1}};
    int fails = 0, load = 1;
    for (auto& c : cases) {
        std::vector<double> T = geom(c.t0, 5, c.S), q(c.S, c.q), L(c.S, c.lam);
        if (c.ramp) for (int t = 0; t < c.S; ++t) if (t >= 0.7 * c.S) L[t] = c.lam * (c.S - 1 - t) / std::max(1.0, c.S - 1 - 0.7 * c.S);
        Tables tb = make_tables(T, q, L, c.jv);
        std::vector<uint64_t> out(size_t(c.trials) * 36);
        std::vector<int32_t> trace(size_t(c.trials) * c.S);
        const uint64_t seed = 0xC0FFEEull; const uint32_t off = 1000;
        sca_v80(J64.data(), tb.fourT.data(), tb.q.data(), tb.kconst.data(), tb.kcorr.data(), out.data(), trace.data(), c.S,
                c.trials, seed, off, g.sumw, load, 1);
        load = 0;  // later cases reuse the BRAM-resident couplings, as on hardware
        for (int tr = 0; tr < c.trials; ++tr) {
            Result r = run_trial(g, tb, seed, off + tr, true);
            bool ok = true;
            for (int w = 0; w < WORDS; ++w) ok &= out[size_t(tr) * 36 + w] == r.spins[w];
            ok &= int64_t(out[size_t(tr) * 36 + 32]) == r.cut && int64_t(out[size_t(tr) * 36 + 33]) == r.flips;
            for (int t = 0; t < c.S; ++t) ok &= trace[size_t(tr) * c.S + t] == r.n_lin[t];
            std::printf("%-20s trial %d: %s (cut %lld, flips %lld)\n", c.name, tr, ok ? "MATCH" : "MISMATCH", (long long)r.cut, (long long)r.flips);
            fails += !ok;
        }
    }
    std::printf(fails ? "TESTBENCH FAILED (%d)\n" : "TESTBENCH PASSED\n", fails);
    return fails ? 1 : 0;
}
