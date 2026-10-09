// V80 Onsager-SCA engine (HLS), v3 (v2 + log-depth priority encoders, per-lane accumulators). Plain C++17 + HLS pragmas, so it also compiles natively for bit-exact tests.
// Semantics are defined by sca_ref.hpp; see ../DESIGN.md.
//
// v2 changes for timing (v1 UPDATE was 8.5 ns) and area (v1: 548k LUT):
//  * One LANE loop (init spins / decide / score) and one ACCUM loop (init fields / update) are shared by all phases.
//  * ACCUM uses two independent flip extractors (spins 0-1023 on BRAM port A, 1024-2047 on port B). Each extracts one
//    bit per cycle and prefetches its next non-empty word, so the loop-carried path is a 64-bit clear-lowest-bit + select.
//  * No diagonal comparators: every spin's own row contributes -d, which keeps the invariant h_stored = h_true - s.
//    Decisions use s*h_true = s*h_stored + 1; the score adds N.
#include <cstdint>

namespace k {
constexpr int N = 2000, NP = 2048, WORDS = NP / 64, HALF = WORDS / 2, LANES = 128, ROUNDS = NP / LANES, MAX_S = 4096;
constexpr int OUT_WORDS = 36;
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
    return v & (~v + 1);  // isolate the lowest set bit: one carry chain
}
static inline int enc64(uint64_t onehot) {  // one-hot -> index with six OR-trees (log depth)
#pragma HLS INLINE
    const uint64_t m[6] = {0xAAAAAAAAAAAAAAAAull, 0xCCCCCCCCCCCCCCCCull, 0xF0F0F0F0F0F0F0F0ull,
                           0xFF00FF00FF00FF00ull, 0xFFFF0000FFFF0000ull, 0xFFFFFFFF00000000ull};
    int n = 0;
    for (int k = 0; k < 6; ++k) {
#pragma HLS UNROLL
        if (onehot & m[k]) n |= 1 << k;
    }
    return n;
}
static inline int enc16(uint32_t onehot) {
#pragma HLS INLINE
    const uint32_t m[4] = {0xAAAAu, 0xCCCCu, 0xF0F0u, 0xFF00u};
    int n = 0;
    for (int k = 0; k < 4; ++k) {
#pragma HLS UNROLL
        if (onehot & m[k]) n |= 1 << k;
    }
    return n;
}
// One flip extractor over 16 words. cur holds the word being drained; occ marks the remaining non-empty words.
struct Extractor {
    uint64_t cur; uint32_t occ; int w;
};
static inline void ex_start(Extractor& e, const uint64_t f[], int base) {
#pragma HLS INLINE
    e.occ = 0;
    for (int j = 0; j < HALF; ++j) {
#pragma HLS UNROLL
        if (f[base + j]) e.occ |= 1u << j;
    }
    e.cur = 0; e.w = 0;
    if (e.occ) {
        const uint32_t lb = e.occ & (~e.occ + 1);
        e.w = enc16(lb); e.cur = f[base + e.w]; e.occ &= ~lb;
    }
}
// Emits one flipped spin per call (if any): index x and its current spin value (from s, after commit).
static inline bool ex_step(Extractor& e, const uint64_t f[], const uint64_t s[], int base, int& x, bool& spin) {
#pragma HLS INLINE
    if (!e.cur) return false;
    const uint64_t lb = lowbit64(e.cur);
    x = (base + e.w) * 64 + enc64(lb);
    spin = (s[base + e.w] & lb) != 0;
    const uint64_t rest = e.cur & ~lb;
    // Prefetch path depends only on occ, not on cur: computed in parallel with the clear-lowest-bit.
    const uint32_t nlb = e.occ & (~e.occ + 1);
    const int nw = enc16(nlb);
    const uint64_t nword = f[base + nw];
    if (rest) {
        e.cur = rest;
    } else if (e.occ) {
        e.cur = nword; e.w = nw; e.occ &= ~nlb;
    } else {
        e.cur = 0;
    }
    return true;
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
    static uint64_t coupling[NP][WORDS];
#pragma HLS ARRAY_PARTITION variable=coupling complete dim=2
#pragma HLS BIND_STORAGE variable=coupling type=ram_t2p impl=bram
    static int32_t fourT[MAX_S], qv[MAX_S], kconst[MAX_S];
    static int64_t kcorr[MAX_S];

    if (load_J) {
    LOAD_J:
        for (int idx = 0; idx < NP * WORDS; ++idx) {
#pragma HLS PIPELINE II=1
            coupling[idx / WORDS][idx % WORDS] = J64[idx];
        }
    }
LOAD_TABLES:
    for (int t = 0; t < S; ++t) {
#pragma HLS PIPELINE II=1
        fourT[t] = fourT_in[t]; qv[t] = q_in[t]; kconst[t] = kconst_in[t]; kcorr[t] = kcorr_in[t];
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
                uint64_t blo = 0, bhi = 0;
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
                    if (l < 64) blo |= b << l; else bhi |= b << (l - 64);
                }
                if (mode == L_INIT) { s[2 * r] = blo; s[2 * r + 1] = bhi; }
                else if (mode == L_DECIDE) { f[2 * r] = blo; f[2 * r + 1] = bhi; }
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
            Extractor ea, eb;
            ex_start(ea, f, 0);
            ex_start(eb, f, HALF);
        ACCUM:
            while (ea.cur || eb.cur) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=1024
                int xa = 0, xb = HALF * 64;
                bool sa = false, sb = false;
                const bool va = ex_step(ea, f, s, 0, xa, sa);
                const bool vb = ex_step(eb, f, s, HALF, xb, sb);
                const int da = va ? (init ? (sa ? 1 : -1) : (sa ? 2 : -2)) : 0;
                const int db = vb ? (init ? (sb ? 1 : -1) : (sb ? 2 : -2)) : 0;
                for (int y = 0; y < NP; ++y) {
#pragma HLS UNROLL
                    const int ca = ((coupling[xa][y >> 6] >> (y & 63)) & 1) ? da : -da;
                    const int cb = ((coupling[xb][y >> 6] >> (y & 63)) & 1) ? db : -db;
                    h[y] += int16_t(ca + cb);
                }
                nflip += uint32_t(va) + uint32_t(vb);
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
