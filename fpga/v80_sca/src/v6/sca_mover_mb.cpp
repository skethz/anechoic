// Multi-bit front end (HLS) for sca_core_mb.v: copy of sca_mover.cpp (v6.2+) with MB_K coupling planes in the J stream.
// Same name, arguments and AXI-Lite register map as sca_mover.cpp, so prepare_ip.py, the timers and the host register
// access are unchanged. J64 holds MB_K planes of 2048 rows x 2048 bits, plane-major (MB_K x 16384 256-bit words); each
// word is streamed as two 128-bit beats, low half first, exactly as before.
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

#ifndef MB_K
#define MB_K 2
#endif

typedef ap_uint<128> beat_t;
typedef ap_uint<64> word_t;

extern "C" void sca_v80(const ap_uint<256>* J64, const ap_uint<256>* tab_in, uint64_t* out, int32_t* trace, int S,
                        int trials, uint64_t seed, uint32_t trial_offset, int64_t sumw, int load_J, int do_trace,
                        hls::stream<beat_t>& to_core, hls::stream<word_t>& from_core) {
#if MB_K == 2
#pragma HLS INTERFACE m_axi port=J64 bundle=gmem0 depth=32768 offset=slave max_read_burst_length=64 num_read_outstanding=16
#elif MB_K == 4
#pragma HLS INTERFACE m_axi port=J64 bundle=gmem0 depth=65536 offset=slave max_read_burst_length=64 num_read_outstanding=16
#elif MB_K == 8
#pragma HLS INTERFACE m_axi port=J64 bundle=gmem0 depth=131072 offset=slave max_read_burst_length=64 num_read_outstanding=16
#else
#error "MB_K must be 2, 4 or 8"
#endif
#pragma HLS INTERFACE m_axi port=tab_in bundle=gmem0 depth=4096 offset=slave max_read_burst_length=64 num_read_outstanding=16
#pragma HLS INTERFACE m_axi port=out bundle=gmem1 depth=36864 offset=slave
#pragma HLS INTERFACE m_axi port=trace bundle=gmem1 depth=65536 offset=slave
#pragma HLS INTERFACE axis port=to_core
#pragma HLS INTERFACE axis port=from_core
#pragma HLS INTERFACE s_axilite port=J64 bundle=control
#pragma HLS INTERFACE s_axilite port=tab_in bundle=control
#pragma HLS INTERFACE s_axilite port=out bundle=control
#pragma HLS INTERFACE s_axilite port=trace bundle=control
#pragma HLS INTERFACE s_axilite port=S bundle=control
#pragma HLS INTERFACE s_axilite port=trials bundle=control
#pragma HLS INTERFACE s_axilite port=seed bundle=control
#pragma HLS INTERFACE s_axilite port=trial_offset bundle=control
#pragma HLS INTERFACE s_axilite port=sumw bundle=control
#pragma HLS INTERFACE s_axilite port=load_J bundle=control
#pragma HLS INTERFACE s_axilite port=do_trace bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
    beat_t h0 = 0, h1 = 0;
    h0.range(15, 0) = ap_uint<16>(trials);
    h0.range(31, 16) = ap_uint<16>(S);
    h0[32] = (load_J & 1) != 0;
    h0[34] = (load_J & 2) != 0;
    h0[33] = do_trace != 0;
    h0.range(127, 64) = seed;
    h1.range(31, 0) = trial_offset;
    h1.range(127, 64) = uint64_t(sumw);
    to_core.write(h0);
    to_core.write(h1);
    if (load_J & 1) {
    LOAD_J:
        for (int i = 0; i < 2048 * 8 * MB_K; ++i) {   // MB_K planes x 2048 rows x 2048 bits -> MB_K x 32768 beats
#pragma HLS PIPELINE II=2
            const ap_uint<256> w = J64[i];
            to_core.write(w.range(127, 0));
            to_core.write(w.range(255, 128));
        }
    }
    if (!(load_J & 2)) {
LOAD_TABLES:
    for (int t = 0; t < S; ++t) {   // {fourT|q, kconst, kcorr, 0} -> {kcorr[31:0], kconst[31:0], q, fourT}
#pragma HLS PIPELINE II=1
        const ap_uint<256> w = tab_in[t];
        beat_t b;
        b.range(63, 0) = w.range(63, 0);
        b.range(95, 64) = w.range(95, 64);
        b.range(127, 96) = w.range(159, 128);
        to_core.write(b);
    }
    }
TRIALS:
    for (int tr = 0; tr < trials; ++tr) {
    RESULT:
        for (int w = 0; w < 36; ++w) {
#pragma HLS PIPELINE II=1
            out[tr * 36 + w] = w < 34 ? uint64_t(from_core.read()) : uint64_t(0);
        }
        if (do_trace) {
        TRACE:
            for (int t = 0; t < S; ++t) {
#pragma HLS PIPELINE II=1
                trace[tr * S + t] = int32_t(uint32_t(from_core.read()));
            }
        }
    }
}
