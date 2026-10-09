// V80 Onsager-SCA engine (HLS), v5.1 (v5 + queue-based extractor refill). Plain C++17 + HLS pragmas; compiles natively.
// Semantics are defined by sca_ref.hpp (build with the same -DSCA_LANES); see ../DESIGN.md. v4 is kept in sca_kernel_v4.cpp.
//
// v5 changes (measured v4 cost: 0.522 cycles/flip + 40.2 cycles/step, and ~35 cycles/step of table loading per launch):
//  * LANES decision lanes per round (v5: 256 -> 8 rounds per step instead of 16).
//  * NEX = 4 flip extractors instead of 2. Extractor e drains words e, e+4, e+8, ... (interleaved for balance) and has its
//    own coupling read port: two identical coupling copies, each true dual-port (extractors 0,1 -> copy 0; 2,3 -> copy 1).
//  * Parameter tables are loaded with four sequential (burstable) loops instead of one interleaved loop.
// Unchanged from v2-v4: one shared LANE loop and one ACCUM loop; diagonal invariant h_stored = h_true - s (decisions use
// s*h_stored + 1, the score adds N); lowest-set-bit extraction with OR-tree encoders and a prefetched next word.
#include <cstdint>

#ifndef SCA_LANES
#define SCA_LANES 128
#endif

namespace k {
constexpr int N = 2000, NP = 2048, WORDS = NP / 64, LANES = SCA_LANES, ROUNDS = NP / LANES, WPR = LANES / 64, MAX_S = 4096;
constexpr int NEX = 4, EXW = WORDS / NEX;  // extractors and words per extractor
constexpr int OUT_WORDS = 36;
static_assert(LANES % 64 == 0 && NP % LANES == 0, "lanes must fill whole words");
static_assert((EXW & (EXW - 1)) == 0, "EXW must be a power of two (queue index mask)");
enum LaneMode { L_INIT = 0, L_DECIDE = 1, L_SCORE = 2 };

struct Lane { uint32_t a, b, c, d; };

static inline uint32_t rotl32(uint32_t x, int r) {
#pragma HLS INLINE
    return (x << r) | (x >> (32 - r));
}
static inline uint32_t xo_next(Lane& L) {
#pragma HLS INLINE
    const uint32_t res = rotl32(L.b * 5u, 7) * 9u, t = L.b << 9;
    L.c ^= L.a; L.d ^= L.b; L.b ^= L.c; L.a ^= L.d; L.c ^= t; L.d = rotl32(L.d, 11);
    return res;
}
static inline uint64_t splitmix(uint64_t& x) {
#pragma HLS INLINE
    uint64_t z = (x += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
static inline uint64_t lowbit64(uint64_t v) {
#pragma HLS INLINE
    return v & (~v + 1);
}
static inline int enc64(uint64_t oh) {  // one-hot -> index; each bit is an independent OR-reduction
#pragma HLS INLINE
    const int b0 = (oh & 0xAAAAAAAAAAAAAAAAull) != 0, b1 = (oh & 0xCCCCCCCCCCCCCCCCull) != 0;
    const int b2 = (oh & 0xF0F0F0F0F0F0F0F0ull) != 0, b3 = (oh & 0xFF00FF00FF00FF00ull) != 0;
    const int b4 = (oh & 0xFFFF0000FFFF0000ull) != 0, b5 = (oh & 0xFFFFFFFF00000000ull) != 0;
    return b0 | (b1 << 1) | (b2 << 2) | (b3 << 3) | (b4 << 4) | (b5 << 5);
}
static inline int enc16(uint32_t oh) {
#pragma HLS INLINE
    const int b0 = (oh & 0xAAAAu) != 0, b1 = (oh & 0xCCCCu) != 0, b2 = (oh & 0xF0F0u) != 0, b3 = (oh & 0xFF00u) != 0;
    return b0 | (b1 << 1) | (b2 << 2) | (b3 << 3);
}
static inline uint32_t lowbit16(uint32_t v) {
#pragma HLS INLINE
    return v & (~v + 1);
}
// Flip extractor over the EXW words {E, E+NEX, E+2*NEX, ...}. v5.1: ex_start compacts the non-empty words into an
// ordered queue (outside the pipelined loop), so a refill is one indexed read with no priority encoder in the loop; the
// step always computes x and spin (no early return), which removes select chains on the coupling address.
struct Extractor {
    int n, head; uint64_t cur; int w;
};
template <int E>
static inline void ex_start(Extractor& e, uint64_t wq[EXW], int wi[EXW], const uint64_t f[]) {
#pragma HLS INLINE
    int n = 0;
    for (int j = 0; j < EXW; ++j) {
#pragma HLS UNROLL
        wq[j] = 0; wi[j] = 0;
    }
    for (int j = 0; j < EXW; ++j) {
#pragma HLS UNROLL
        const uint64_t v = f[j * NEX + E];
        if (v) { wq[n] = v; wi[n] = j; ++n; }
    }
    e.n = n; e.head = 1; e.cur = wq[0]; e.w = wi[0];
}
template <int E>
static inline bool ex_step(Extractor& e, const uint64_t wq[EXW], const int wi[EXW], const uint64_t s[], int& x, bool& spin) {
#pragma HLS INLINE
    const bool valid = e.cur != 0;
    const uint64_t lb = lowbit64(e.cur);
    const int word = e.w * NEX + E;
    x = (word << 6) | enc64(lb);
    spin = (s[word] & lb) != 0;
    const uint64_t rest = e.cur & ~lb;
    const uint64_t qw = wq[e.head & (EXW - 1)];
    const int qi = wi[e.head & (EXW - 1)];
    const bool more = e.head < e.n;
    if (rest) {
        e.cur = rest;
    } else {
        e.cur = more ? qw : 0;
        e.w = qi;
        e.head = e.head + 1;
    }
    return valid;
}
}  // namespace k

extern "C" void sca_v80(const uint64_t* J64, const int32_t* fourT_in, const int32_t* q_in, const int32_t* kconst_in,
                        const int64_t* kcorr_in, uint64_t* out, int32_t* trace, int S, int trials, uint64_t seed,
                        uint32_t trial_offset, int64_t sumw, int load_J, int do_trace) {
#pragma HLS INTERFACE m_axi port=J64 bundle=gmem0 depth=65536 offset=slave
#pragma HLS INTERFACE m_axi port=fourT_in bundle=gmem1 depth=4096 offset=slave
#pragma HLS INTERFACE m_axi port=q_in bundle=gmem1 depth=4096 offset=slave
#pragma HLS INTERFACE m_axi port=kconst_in bundle=gmem1 depth=4096 offset=slave
#pragma HLS INTERFACE m_axi port=kcorr_in bundle=gmem1 depth=4096 offset=slave
#pragma HLS INTERFACE m_axi port=out bundle=gmem2 depth=36864 offset=slave
#pragma HLS INTERFACE m_axi port=trace bundle=gmem2 depth=65536 offset=slave
#pragma HLS INTERFACE s_axilite port=J64 bundle=control
#pragma HLS INTERFACE s_axilite port=fourT_in bundle=control
#pragma HLS INTERFACE s_axilite port=q_in bundle=control
#pragma HLS INTERFACE s_axilite port=kconst_in bundle=control
#pragma HLS INTERFACE s_axilite port=kcorr_in bundle=control
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
    using namespace k;
    static uint64_t coupling0[NP][WORDS], coupling1[NP][WORDS];   // identical copies: 4 read ports in total
#pragma HLS ARRAY_PARTITION variable=coupling0 complete dim=2
#pragma HLS ARRAY_PARTITION variable=coupling1 complete dim=2
#pragma HLS BIND_STORAGE variable=coupling0 type=ram_t2p impl=bram
#pragma HLS BIND_STORAGE variable=coupling1 type=ram_t2p impl=bram
    static int32_t fourT[MAX_S], qv[MAX_S], kconst[MAX_S];
    static int64_t kcorr[MAX_S];

    if (load_J) {
    LOAD_J:
        for (int idx = 0; idx < NP * WORDS; ++idx) {
#pragma HLS PIPELINE II=1
            const uint64_t v = J64[idx];
            coupling0[idx / WORDS][idx % WORDS] = v;
            coupling1[idx / WORDS][idx % WORDS] = v;
        }
    }
LOAD_FOURT:
    for (int t = 0; t < S; ++t) {
#pragma HLS PIPELINE II=1
        fourT[t] = fourT_in[t];
    }
LOAD_Q:
    for (int t = 0; t < S; ++t) {
#pragma HLS PIPELINE II=1
        qv[t] = q_in[t];
    }
LOAD_KCONST:
    for (int t = 0; t < S; ++t) {
#pragma HLS PIPELINE II=1
        kconst[t] = kconst_in[t];
    }
LOAD_KCORR:
    for (int t = 0; t < S; ++t) {
#pragma HLS PIPELINE II=1
        kcorr[t] = kcorr_in[t];
    }

TRIALS:
    for (int tr = 0; tr < trials; ++tr) {
        const uint32_t trial = trial_offset + uint32_t(tr);
        Lane rng[LANES];
#pragma HLS ARRAY_PARTITION variable=rng complete
        uint64_t s[WORDS], sp[WORDS], f[WORDS];
#pragma HLS ARRAY_PARTITION variable=s complete
#pragma HLS ARRAY_PARTITION variable=sp complete
#pragma HLS ARRAY_PARTITION variable=f complete
        int16_t h[NP];  // stored field = true field - s (see header)
#pragma HLS ARRAY_PARTITION variable=h complete

    SEED:
        for (int l = 0; l < LANES; ++l) {
#pragma HLS PIPELINE II=1
            uint64_t x = seed ^ (uint64_t(trial) << 20) ^ (uint64_t(l) << 48) ^ 0x5CA0F00Dull;
            const uint64_t a = splitmix(x), b = splitmix(x);
            Lane L{uint32_t(a), uint32_t(a >> 32), uint32_t(b), uint32_t(b >> 32)};
            if (!(L.a | L.b | L.c | L.d)) L.a = 1;
            rng[l] = L;
        }
        for (int w = 0; w < WORDS; ++w) {
#pragma HLS UNROLL
            sp[w] = ~0ull; s[w] = 0;
        }
        for (int y = 0; y < NP; ++y) {
#pragma HLS UNROLL
            h[y] = 0;
        }
        int32_t nlin_prev = 0;
        int64_t flips_total = 0, shsum = 0;

        // Phase t = -1: init spins (LANE) + init fields (ACCUM); t = 0..S-1: decide + update; t = S: score.
    PHASES:
        for (int t = -1; t <= S; ++t) {
            const int mode = t < 0 ? L_INIT : (t < S ? L_DECIDE : L_SCORE);
            const int tt = t < 0 ? 0 : (t < S ? t : 0);
            const int32_t corr = int32_t(((int64_t(nlin_prev) * kcorr[tt]) >> 8) + kconst[tt]);
            const int32_t ft = fourT[tt], twoT = ft >> 1, qt = qv[tt];
            int32_t lacc[LANES];
#pragma HLS ARRAY_PARTITION variable=lacc complete
            for (int l = 0; l < LANES; ++l) {
#pragma HLS UNROLL
                lacc[l] = 0;
            }
        LANE:
            for (int r = 0; r < ROUNDS; ++r) {
#pragma HLS PIPELINE II=1
                uint64_t bw[WPR];
#pragma HLS ARRAY_PARTITION variable=bw complete
                for (int k2 = 0; k2 < WPR; ++k2) {
#pragma HLS UNROLL
                    bw[k2] = 0;
                }
                for (int l = 0; l < LANES; ++l) {
#pragma HLS UNROLL
                    const int i = r * LANES + l;
                    const bool valid = i < N;
                    const int si = int((s[i >> 6] >> (i & 63)) & 1), spi = int((sp[i >> 6] >> (i & 63)) & 1);
                    const int32_t sh = (si ? int32_t(h[i]) : -int32_t(h[i])) + 1;  // s_i * h_true
                    uint64_t b = 0;
                    if (mode != L_SCORE) {
                        const uint32_t rv = xo_next(rng[l]);
                        if (mode == L_INIT) {
                            b = valid ? uint64_t(rv >> 31) : 1ull;  // padding spins fixed at +1
                        } else {
                            const int32_t z = sh * 65536 + qt - (si == spi ? corr : -corr);
                            const int32_t thr = int32_t((int64_t(int32_t(rv >> 16) - 32768) * ft) >> 16);
                            b = (valid && z < thr) ? 1ull : 0ull;
                            lacc[l] += (valid && z > -twoT && z < twoT) ? 1 : 0;  // per-lane n_lin count
                        }
                    } else if (valid) {
                        lacc[l] += sh;  // per-lane score partial
                    }
                    bw[l >> 6] |= b << (l & 63);
                }
                for (int k2 = 0; k2 < WPR; ++k2) {
#pragma HLS UNROLL
                    if (mode == L_INIT) s[r * WPR + k2] = bw[k2];
                    else if (mode == L_DECIDE) f[r * WPR + k2] = bw[k2];
                }
            }
            int64_t lsum = 0;  // reduce the per-lane accumulators once per phase (outside the II=1 loop)
            for (int l = 0; l < LANES; ++l) {
#pragma HLS UNROLL
                lsum += lacc[l];
            }
            const int32_t nlin = int32_t(lsum);
            if (mode == L_SCORE) { shsum = lsum; break; }

            // ACCUM: init mode adds d = s_x for every active spin; update mode adds d = 2*s_x(new) for flipped spins.
            const bool init = mode == L_INIT;
            uint32_t nflip = 0;
            for (int w = 0; w < WORDS; ++w) {
#pragma HLS UNROLL
                if (init) {
                    f[w] = (w == WORDS - 1) ? ((1ull << (N - 64 * (WORDS - 1))) - 1) : ~0ull;  // all active spins
                } else {
                    sp[w] = s[w];
                    s[w] ^= f[w];
                }
            }
            Extractor e0, e1, e2, e3;
            uint64_t q0[EXW], q1[EXW], q2[EXW], q3[EXW];
            int i0[EXW], i1[EXW], i2[EXW], i3[EXW];
#pragma HLS ARRAY_PARTITION variable=q0 complete
#pragma HLS ARRAY_PARTITION variable=q1 complete
#pragma HLS ARRAY_PARTITION variable=q2 complete
#pragma HLS ARRAY_PARTITION variable=q3 complete
#pragma HLS ARRAY_PARTITION variable=i0 complete
#pragma HLS ARRAY_PARTITION variable=i1 complete
#pragma HLS ARRAY_PARTITION variable=i2 complete
#pragma HLS ARRAY_PARTITION variable=i3 complete
            ex_start<0>(e0, q0, i0, f); ex_start<1>(e1, q1, i1, f); ex_start<2>(e2, q2, i2, f); ex_start<3>(e3, q3, i3, f);
        ACCUM:
            while (e0.cur || e1.cur || e2.cur || e3.cur) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=512
                int x0 = 0, x1 = 1 << 6, x2 = 2 << 6, x3 = 3 << 6;
                bool s0 = false, s1 = false, s2 = false, s3 = false;
                const bool v0 = ex_step<0>(e0, q0, i0, s, x0, s0);
                const bool v1 = ex_step<1>(e1, q1, i1, s, x1, s1);
                const bool v2 = ex_step<2>(e2, q2, i2, s, x2, s2);
                const bool v3 = ex_step<3>(e3, q3, i3, s, x3, s3);
                const int m = init ? 1 : 2;
                const int d0 = v0 ? (s0 ? m : -m) : 0, d1 = v1 ? (s1 ? m : -m) : 0;
                const int d2 = v2 ? (s2 ? m : -m) : 0, d3 = v3 ? (s3 ? m : -m) : 0;
                for (int y = 0; y < NP; ++y) {
#pragma HLS UNROLL
                    const int c0 = ((coupling0[x0][y >> 6] >> (y & 63)) & 1) ? d0 : -d0;
                    const int c1 = ((coupling0[x1][y >> 6] >> (y & 63)) & 1) ? d1 : -d1;
                    const int c2 = ((coupling1[x2][y >> 6] >> (y & 63)) & 1) ? d2 : -d2;
                    const int c3 = ((coupling1[x3][y >> 6] >> (y & 63)) & 1) ? d3 : -d3;
                    h[y] += int16_t((c0 + c1) + (c2 + c3));
                }
                nflip += uint32_t(v0) + uint32_t(v1) + uint32_t(v2) + uint32_t(v3);
            }
            if (!init) {
                flips_total += nflip;
                nlin_prev = nlin;
                if (do_trace) trace[int64_t(tr) * S + t] = nlin;
            }
        }
        const int64_t base = int64_t(tr) * OUT_WORDS;
    WRITE:
        for (int w = 0; w < WORDS; ++w) {
#pragma HLS PIPELINE II=1
            uint64_t v = s[w];
            if (w == WORDS - 1) v &= (1ull << (N - 64 * (WORDS - 1))) - 1;
            out[base + w] = v;
        }
        out[base + 32] = uint64_t((sumw + shsum / 2) / 2);
        out[base + 33] = uint64_t(flips_total);
        out[base + 34] = 0;
        out[base + 35] = 0;
    }
}
