// v6 kernel front end (HLS): same name, arguments and AXI-Lite register map as sca_v80 v5.x, so the host, timers and
// protocol scripts are unchanged. It only moves data: header, couplings (optional) and tables into the RTL core
// (sca_core.v, 128-bit AXI-Stream), then the per-trial results back to memory. All computation is in the core.
// Memory: J64 and tab_in through one 256-bit read master (gmem0); out and trace through one 64-bit write master (gmem1),
// so each engine uses two NoC/SmartConnect ports instead of three.
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

typedef ap_uint<128> beat_t;
typedef ap_uint<64> word_t;

extern "C" void sca_v80(const ap_uint<256>* J64, const ap_uint<256>* tab_in, uint64_t* out, int32_t* trace, int S,
                        int trials, uint64_t seed, uint32_t trial_offset, int64_t sumw, int load_J, int do_trace,
                        hls::stream<beat_t>& to_core, hls::stream<word_t>& from_core) {
#pragma HLS INTERFACE m_axi port=J64 bundle=gmem0 depth=16384 offset=slave max_read_burst_length=64 num_read_outstanding=16
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
    h0[32] = load_J != 0;
    h0[33] = do_trace != 0;
    h0.range(127, 64) = seed;
    h1.range(31, 0) = trial_offset;
    h1.range(127, 64) = uint64_t(sumw);
    to_core.write(h0);
    to_core.write(h1);
    if (load_J) {
    LOAD_J:
        for (int i = 0; i < 2048 * 8; ++i) {   // 2048 rows x 2048 bits = 16384 x 256-bit words -> 32768 beats
#pragma HLS PIPELINE II=2
            const ap_uint<256> w = J64[i];
            to_core.write(w.range(127, 0));
            to_core.write(w.range(255, 128));
        }
    }
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
