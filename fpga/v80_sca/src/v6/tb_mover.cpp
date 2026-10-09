// Unit testbench for the v6 front end (sca_mover.cpp): checks the exact beat sequence sent to the core and the memory
// image written from the core's result stream. The core itself is verified separately (tb_core.v, bit-exact vs
// sca_ref.hpp); the assembled engine is verified on the board by the protocol's exact-reference gate.
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern "C" void sca_v80(const ap_uint<256>*, const ap_uint<256>*, uint64_t*, int32_t*, int, int, uint64_t, uint32_t,
                        int64_t, int, int, hls::stream<ap_uint<128>>&, hls::stream<ap_uint<64>>&);

static uint64_t rng_state = 0x9E3779B97F4A7C15ull;
static uint64_t rnd() { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

int main() {
    std::vector<ap_uint<256>> J(16384), tab(4096);
    std::vector<uint64_t> Jw(65536);
    for (int i = 0; i < 65536; ++i) Jw[i] = rnd();
    for (int i = 0; i < 16384; ++i)
        for (int k = 0; k < 4; ++k) J[i].range(64 * k + 63, 64 * k) = Jw[4 * i + k];
    std::vector<uint64_t> out(36864, 0xDEADull);
    std::vector<int32_t> trace(65536, -7);
    int fails = 0;
    struct Case { int S, trials, loadJ, trace; };   // loadJ bit 0: couplings, bit 1: keep tables (v6.2)
    const Case cases[] = {{37, 3, 1, 1}, {5, 2, 0, 0}, {5, 2, 2, 1}, {9, 1, 3, 0}};
    for (const Case& c : cases) {
        std::vector<uint64_t> tw(4 * c.S);
        for (int t = 0; t < c.S; ++t) {
            for (int k = 0; k < 4; ++k) tw[4 * t + k] = rnd();
            for (int k = 0; k < 4; ++k) tab[t].range(64 * k + 63, 64 * k) = tw[4 * t + k];
        }
        hls::stream<ap_uint<128>> to_core("to_core");
        hls::stream<ap_uint<64>> from_core("from_core");
        const int per = 34 + (c.trace ? c.S : 0);
        std::vector<uint64_t> res(size_t(c.trials) * per);
        for (auto& v : res) { v = rnd(); from_core.write(v); }
        const uint64_t seed = rnd(); const uint32_t off = uint32_t(rnd()); const int64_t sumw = int64_t(rnd());
        sca_v80(J.data(), tab.data(), out.data(), trace.data(), c.S, c.trials, seed, off, sumw, c.loadJ, c.trace, to_core,
                from_core);
        // header
        ap_uint<128> h0 = to_core.read(), h1 = to_core.read();
        fails += h0.range(15, 0) != unsigned(c.trials) || h0.range(31, 16) != unsigned(c.S) || h0[32] != (c.loadJ & 1) ||
                 h0[34] != ((c.loadJ >> 1) & 1) ||
                 h0[33] != c.trace || uint64_t(h0.range(127, 64)) != seed;
        fails += uint32_t(h1.range(31, 0)) != off || uint64_t(h1.range(127, 64)) != uint64_t(sumw);
        if (c.loadJ & 1)
            for (int b = 0; b < 32768; ++b) {
                ap_uint<128> v = to_core.read();
                fails += uint64_t(v.range(63, 0)) != Jw[2 * b] || uint64_t(v.range(127, 64)) != Jw[2 * b + 1];
            }
        for (int t = 0; t < ((c.loadJ & 2) ? 0 : c.S); ++t) {
            ap_uint<128> v = to_core.read();
            fails += uint32_t(v.range(31, 0)) != uint32_t(tw[4 * t]) || uint32_t(v.range(63, 32)) != uint32_t(tw[4 * t] >> 32);
            fails += uint32_t(v.range(95, 64)) != uint32_t(tw[4 * t + 1]) || uint32_t(v.range(127, 96)) != uint32_t(tw[4 * t + 2]);
        }
        fails += !to_core.empty() || !from_core.empty();
        for (int tr = 0; tr < c.trials; ++tr) {
            for (int w = 0; w < 36; ++w) fails += out[tr * 36 + w] != (w < 34 ? res[size_t(tr) * per + w] : 0ull);
            if (c.trace)
                for (int t = 0; t < c.S; ++t) fails += trace[tr * c.S + t] != int32_t(uint32_t(res[size_t(tr) * per + 34 + t]));
        }
        std::printf("case S=%d trials=%d loadJ=%d trace=%d: %s\n", c.S, c.trials, c.loadJ, c.trace, fails ? "FAIL" : "ok");
    }
    std::printf(fails ? "TESTBENCH FAILED (%d)\n" : "TESTBENCH PASSED\n", fails);
    return fails ? 1 : 0;
}
