// K-bit coupling + bias (external field) version of the bit-exact GH200 SCA engine.
// Reference: sca::run_trial_bias(n, dense, bias, ...) (fpga/v80_sca/src/sca_ref_bias.hpp, SCA_LANES=256); arithmetic:
// fpga/v80_sca/MULTIBIT_SPEC.md incl. the section "Bias (external field)". Derived from gpu/multibit_20261007/src/sca_gpu_mb.cu
// (left untouched); the only kernel changes are
//   (1) the fields start at the bias, h(0) = b + J s(0) (the initial-field BMMA pass is unchanged; the per-step update is
//       unchanged, so the bias persists exactly), and
//   (2) an optional WIDE instantiation that clamps |h| to hclamp inside the decision. It is used only when the host cannot
//       prove that s*h*65536 + q -/+ corr fits int32 (large biases); clamping is exact (host::plan_ranges).
// The number of active spins n (1..2048) is a runtime parameter (slots >= n are padding: never decided, b = J = 0, s = +1;
// every RNG lane still consumes all 8 draws per step, as in the reference).
//
// Couplings J_xy in {-M..M}, M = 2^(K-1) - 1, stored as K bit planes: plane 0 = sign (1 iff J < 0), plane 1+b = bit b of |J|.
// For the flips F of a step (new spins sigma'), Delta h_y = 2 sum_{x in F} sigma'_x J_xy. With binary tensor-core AND+POPC
// (mma.m16n8k256 .b1 .and.popc), two BMMAs per magnitude plane b give, per field slot y and chain,
//     cA1 = popc(m_b & F+),  cB1 = popc(m_b & F),  cA2 = popc(m_b & sg & F+),  cB2 = popc(m_b & sg & F)     (F+ = F & [sigma' = +1])
// and the spec's two population counts are  popc_e(m_b & [tau_e = +1]) = cA1 + cB2 - 2 cA2  and  popc_e(m_b) = cB1, so
//     sum_{x in F} sigma'_x J_xy = sum_b 2^b (2 cA1 - cB1 - 4 cA2 + 2 cB2)          (exact integer identity; J_yy = 0, padding 0).
// m_b & sg is formed from the resident planes with LOP3 before each BMMA (no extra storage).
//
// Layout (template R): the cluster has 32R warps, warp wg owns 8/R RNG lanes x 8 rounds = 64/R field slots and the full x range.
//   R = 1: 32 warps (v2 mapping, 8 slots per lane). R = 2: 64 warps, 4 slots per lane (two lanes share an RNG lane stream, each
//   generates all 8 draws and uses its 4 rounds). Planes 0..PREG-1 are register-resident A fragments, planes PREG..K-1 live in
//   shared memory (read every step). K = 2 with R = 2, PREG = 2 keeps J (1 MB) in registers of 8 CTAs x 8 warps (one CTA per SM).
// Per step and group: decisions (exactly as run_trial, Q16.16 integer) -> per-chain masks F, F+ (and n_lin counts when needed)
// -> st.async (DSMEM, mbarrier complete_tx) to every CTA -> BMMAs -> exact field update h += 2 * sum (factor 1 for the initial field).
#include "sca_ref.hpp"
#include "tables.hpp"
#include "jmat.hpp"
#include <cuda_runtime.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if SCA_LANES != 256
#error "build with -DSCA_LANES=256"
#endif

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) throw std::runtime_error(std::string(#x) + ": " + cudaGetErrorString(e_)); } while (0)

static const char* kExpectedUuid = "8f88f1cd-2cd1-ab46-b984-2e97ebd18861";

namespace kern {

constexpr int NPAD = 2048;

struct Params {
    const uint4* A;          // [32R warps][K planes][JB][8 kb][32 lanes] uint4
    const int32_t* fourT;    // [S]
    const int32_t* q;        // [S]
    const int64_t* kcorr;    // [S]
    const int32_t* kconst;   // [S]
    int S;
    int N;                   // active spins (slots >= N are padding: never decided, couplings 0)
    const int32_t* bias;     // [2048] per-slot bias b (0 on padding), or nullptr for b = 0
    int32_t hclamp;          // WIDE kernels only: |h| is clamped to hclamp in the decision (exact; host::plan_ranges)
    uint64_t seed;
    uint32_t trial_offset;
    int nchains;             // chains reported (chain index < nchains); other chains of a group are computed but not written
    int do_trace;
    uint8_t* out_spins;      // [nchains][2048]: 1 iff the slot's final spin is +1
    int32_t* out_h;          // [nchains][2048]: final field h of the slot
    int32_t* out_flips;      // [nchains][32R][8]: per-lane partial flip counts
    int32_t* out_trace;      // [nchains][S] n_lin per step
    long long* prof;         // [gridDim][8] phase cycle totals (SCA_PROFILE builds only)
};

__device__ __forceinline__ uint32_t smem_u32(const void* p) { return (uint32_t)__cvta_generic_to_shared(p); }
__device__ __forceinline__ uint32_t mapa(uint32_t a, uint32_t rank) {
    uint32_t r; asm volatile("mapa.shared::cluster.u32 %0, %1, %2;" : "=r"(r) : "r"(a), "r"(rank)); return r;
}
__device__ __forceinline__ uint32_t cluster_rank() { uint32_t r; asm volatile("mov.u32 %0, %%cluster_ctarank;" : "=r"(r)); return r; }
__device__ __forceinline__ uint32_t cluster_index() { uint32_t r; asm volatile("mov.u32 %0, %%clusterid.x;" : "=r"(r)); return r; }
__device__ __forceinline__ void cluster_sync_all() {
    asm volatile("barrier.cluster.arrive.release.aligned;\n\tbarrier.cluster.wait.acquire.aligned;" ::: "memory");
}
__device__ __forceinline__ void mbar_init(uint32_t bar, uint32_t count) {
    asm volatile("mbarrier.init.shared::cta.b64 [%0], %1;" :: "r"(bar), "r"(count) : "memory");
}
__device__ __forceinline__ void mbar_arm(uint32_t bar, uint32_t bytes) {
    asm volatile("mbarrier.arrive.expect_tx.release.cta.shared::cta.b64 _, [%0], %1;" :: "r"(bar), "r"(bytes) : "memory");
}
__device__ __forceinline__ bool mbar_test(uint32_t bar, uint32_t parity) {
    uint32_t ok;
    asm volatile("{\n\t.reg .pred p;\n\tmbarrier.try_wait.parity.acquire.cluster.shared::cta.b64 p, [%1], %2;\n\tselp.u32 %0, 1, 0, p;\n\t}"
                 : "=r"(ok) : "r"(bar), "r"(parity) : "memory");
    return ok != 0;
}
__device__ __forceinline__ void st_async8(uint32_t raddr, uint32_t rbar, uint32_t a, uint32_t b) {
    asm volatile("st.async.shared::cluster.mbarrier::complete_tx::bytes.v2.b32 [%0], {%1,%2}, [%3];"
                 :: "r"(raddr), "r"(a), "r"(b), "r"(rbar) : "memory");
}
__device__ __forceinline__ void st_async16(uint32_t raddr, uint32_t rbar, uint4 v) {
    asm volatile("st.async.shared::cluster.mbarrier::complete_tx::bytes.v4.b32 [%0], {%1,%2,%3,%4}, [%5];"
                 :: "r"(raddr), "r"(v.x), "r"(v.y), "r"(v.z), "r"(v.w), "r"(rbar) : "memory");
}
__device__ __forceinline__ void st_async4(uint32_t raddr, uint32_t rbar, uint32_t v) {
    asm volatile("st.async.shared::cluster.mbarrier::complete_tx::bytes.b32 [%0], %1, [%2];" :: "r"(raddr), "r"(v), "r"(rbar) : "memory");
}
__device__ __forceinline__ void bmma(int32_t (&c)[4], const uint4& a, uint32_t b0, uint32_t b1) {
    asm("mma.sync.aligned.m16n8k256.row.col.s32.b1.b1.s32.and.popc {%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%0,%1,%2,%3};"
        : "+r"(c[0]), "+r"(c[1]), "+r"(c[2]), "+r"(c[3]) : "r"(a.x), "r"(a.y), "r"(a.z), "r"(a.w), "r"(b0), "r"(b1));
}

// xoshiro128**, seeded by splitmix64(seed, trial, lane) -- identical to sca::Rng (copied from the v2 kernel).
struct Rng { uint32_t s0, s1, s2, s3; };
__device__ __forceinline__ uint64_t splitmix(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
__device__ __forceinline__ void rng_seed(Rng& r, uint64_t seed, uint32_t trial, uint32_t lane) {
    uint64_t x = seed ^ (uint64_t(trial) << 20) ^ (uint64_t(lane) << 48) ^ 0x5CA0F00Dull;
    const uint64_t a = splitmix(x), b = splitmix(x);
    r.s0 = uint32_t(a); r.s1 = uint32_t(a >> 32); r.s2 = uint32_t(b); r.s3 = uint32_t(b >> 32);
    if (!(r.s0 | r.s1 | r.s2 | r.s3)) r.s0 = 1;
}
__device__ __forceinline__ uint32_t rng_next(Rng& r) {
    const uint32_t m = r.s1 * 5u;
    const uint32_t res = __funnelshift_l(m, m, 7) * 9u, t = r.s1 << 9;
    r.s2 ^= r.s0; r.s3 ^= r.s1; r.s1 ^= r.s2; r.s0 ^= r.s3; r.s2 ^= t; r.s3 = __funnelshift_l(r.s3, r.s3, 11);
    return res;
}
// bit q of a 4-bit value -> bit 2q
__device__ __forceinline__ uint32_t spread4(uint32_t x) { return (x & 1u) | ((x & 2u) << 1) | ((x & 4u) << 2) | ((x & 8u) << 3); }

template <int K, int R, int PREG, int CS, int G>
struct Cfg {
    static constexpr int NW = 32 * R;          // warps per cluster
    static constexpr int NWG = NW / CS;        // warps per CTA
    static constexpr int THREADS = 32 * NWG;
    static constexpr int LPW = 8 / R;          // RNG lanes per warp
    static constexpr int JB = 4 / R;           // m16 row blocks per warp
    static constexpr int NSL = 8 / R;          // slots (decisions) per lane and chain
    static constexpr int MAGS = K - 1;         // magnitude planes
    static constexpr int PSM = K - PREG;       // planes in shared memory
    static constexpr int KS = (MAGS * 2 * JB < 8 && G <= 2) ? 2 : 1;   // accumulator split over kb parity (more independent BMMA chains)
    // G >= 4: the other groups cover the message flight, so the next draws fill the tensor-pipe shadow instead (measured better
    // at G = 4, worse at G = 2, 3; development record in sca_gpu_mb_v3_dev.cu).
    static constexpr bool RNG_IN_RECV = G >= 4;
    // Message (per group and step): B words [64 W][8 columns n = 2c + kind] (kind 0 = F+, 1 = F), then one n_lin word per
    // source warp (byte c = chain c).
    static constexpr int MASK_WORDS = 64 * 8;
    static constexpr int CNT_WORDS = NW;
    static constexpr int MW = MASK_WORDS + CNT_WORDS;
    static constexpr int RB_OFF = 0;                                        // [G][2][MW] u32
    static constexpr int BAR_OFF = RB_OFF + G * 2 * MW * 4;                 // [G][2] u64
    static constexpr int PL_OFF = ((BAR_OFF + G * 2 * 8) + 15) / 16 * 16;  // [NWG][PSM][JB][8][32] uint4
    static constexpr int PL_BYTES = NWG * PSM * JB * 8 * 32 * 16;
    static constexpr int TAB_OFF = PL_OFF + PL_BYTES;                       // fourT[S], q[S], kconst[S] (i32), kcorr[S] (i64)
    static_assert(R == 1 || R == 2, "R must be 1 or 2");
    static_assert(NW % CS == 0 && CS >= 4 && CS <= 16, "cluster size");
    static_assert(PREG >= 1 && PREG <= K && K >= 2, "planes");
    static_assert(TAB_OFF % 16 == 0 && (MW * 4) % 16 == 0, "alignment");
};

// ONS: the correction depends on the previous step's n_lin. Otherwise corr = kconst[t] is a per-step constant.
// LIN: count and exchange n_lin (needed for ONS and for the trace).
// No CTA-wide barrier inside the step loop (as v2): buffer (group, parity) is rewritten with message m+2 only after every warp
// of the cluster has sent message m+1, i.e. after every warp has finished reading message m.
// Remote stores per (source warp, destination CTA) and step: four 8-byte (R = 2) or four 16-byte (R = 1) mask stores + one
// packed n_lin word (LIN builds). A first version sent the n_lin counts per chain (8 remote stores per warp and CTA) and was
// limited by the receivers' remote-store rate (sca_gpu_mb_v0_dev.cu).
template <int K, int R, int PREG, int CS, int G, bool ONS, bool LIN, bool WIDE>
__global__ void __launch_bounds__(Cfg<K, R, PREG, CS, G>::THREADS, 1) mb_kernel(const Params p) {
    static_assert(LIN || !ONS, "Onsager mode needs n_lin");
    using C = Cfg<K, R, PREG, CS, G>;
    constexpr int NW = C::NW, NWG = C::NWG, LPW = C::LPW, JB = C::JB, NSL = C::NSL, MAGS = C::MAGS, PSM = C::PSM, KS = C::KS;
    constexpr int MW = C::MW, MASKW = C::MASK_WORDS;
    constexpr uint32_t BYTES = uint32_t((C::MASK_WORDS + (LIN ? C::CNT_WORDS : 0)) * 4);
    constexpr uint32_t ALLS = (1u << NSL) - 1u;
    constexpr uint32_t FULL = 0xffffffffu;
    extern __shared__ __align__(16) uint8_t smem[];
    uint32_t* RB = reinterpret_cast<uint32_t*>(smem + C::RB_OFF);
    uint4* PLs = reinterpret_cast<uint4*>(smem + C::PL_OFF);
    const int S = p.S;
    int32_t* T_fourT = reinterpret_cast<int32_t*>(smem + C::TAB_OFF);
    int32_t* T_q = T_fourT + S;
    int32_t* T_kconst = T_q + S;
    int64_t* T_kcorr = reinterpret_cast<int64_t*>(smem + C::TAB_OFF + ((12 * S + 15) / 16) * 16);

    const int tid = threadIdx.x, L = tid & 31, w = tid >> 5;
    const uint32_t rank = cluster_rank();
    const int wg = int(rank) * NWG + w;                 // warp-global index 0..NW-1
    const int c4 = L & 3;                               // chain within group (C-fragment column pair)
    const int lam = (L >> 2) % LPW;                     // RNG lane offset within the warp
    const int hsel = (L >> 2) / LPW;                    // R = 2: which 4 of the 8 rounds (0: even, 1: odd); R = 1: 0
    const int lane_l = LPW * wg + lam;                  // RNG lane 0..255
    // local slot q: round = q (R = 1) or 2q + hsel (R = 2); slot index = round*256 + lane_l
    uint32_t validm = 0;
#pragma unroll
    for (int q = 0; q < NSL; ++q)
        if ((R == 1 ? q : 2 * q + hsel) * 256 + lane_l < p.N) validm |= 1u << q;
    const int cl = int(cluster_index());
#ifdef SCA_PROFILE
    long long pr[8] = {0, 0, 0, 0, 0, 0, 0, 0}; long long tq = clock64();
#define PROF(i) do { const long long tn = clock64(); pr[i] += tn - tq; tq = tn; } while (0)
#else
#define PROF(i) do {} while (0)
#endif

    for (int i = tid; i < S; i += blockDim.x) {
        T_fourT[i] = p.fourT[i]; T_q[i] = p.q[i]; T_kconst[i] = p.kconst[i]; T_kcorr[i] = p.kcorr[i];
    }
    const uint32_t bar0 = smem_u32(smem + C::BAR_OFF);
    if (tid == 0) {
        for (int i = 0; i < 2 * G; ++i) mbar_init(bar0 + 8 * i, 1);
        asm volatile("fence.mbarrier_init.release.cluster;" ::: "memory");
        for (int i = 0; i < 2 * G; ++i) mbar_arm(bar0 + 8 * i, BYTES);   // messages 0 and 1 of every group
    }
    // A fragments: planes 0..PREG-1 resident in registers for the whole kernel, planes PREG..K-1 in shared memory
    uint4 A[PREG][JB][8];
#pragma unroll
    for (int pl = 0; pl < PREG; ++pl)
#pragma unroll
        for (int j = 0; j < JB; ++j)
#pragma unroll
            for (int kb = 0; kb < 8; ++kb) A[pl][j][kb] = p.A[(((wg * K + pl) * JB + j) * 8 + kb) * 32 + L];
    if constexpr (PSM > 0) {
#pragma unroll 1
        for (int pp = 0; pp < PSM; ++pp)
#pragma unroll 1
            for (int j = 0; j < JB; ++j)
#pragma unroll
                for (int kb = 0; kb < 8; ++kb)
                    PLs[(((w * PSM + pp) * JB + j) * 8 + kb) * 32 + L] = p.A[(((wg * K + PREG + pp) * JB + j) * 8 + kb) * 32 + L];
    }
    auto plane = [&](int pl, int j, int kb) -> uint4 {   // pl, j, kb are compile-time after unrolling
        if (pl < PREG) return A[pl < PREG ? pl : 0][j][kb];
        return PLs[(((w * PSM + (pl - PREG)) * JB + j) * 8 + kb) * 32 + L];
    };

    // draws of one step for this lane's slots (all 8 draws of the RNG lane are consumed, as in run_trial)
    auto draw_slots = [&](Rng& r, uint32_t (&d)[NSL]) {
        if constexpr (R == 1) {
#pragma unroll
            for (int k = 0; k < 8; ++k) d[k] = rng_next(r);
        } else {
#pragma unroll
            for (int q = 0; q < 4; ++q) { const uint32_t d0 = rng_next(r), d1 = rng_next(r); d[q] = hsel ? d1 : d0; }
        }
    };

    // bias of this lane's slots: the fields start at h(0) = b + J s(0) (the J s(0) part comes from message 0 below)
    int32_t b0[NSL];
#pragma unroll
    for (int q = 0; q < NSL; ++q) b0[q] = p.bias ? p.bias[(R == 1 ? q : 2 * q + hsel) * 256 + lane_l] : 0;

    // per-group chain state; ntk / l1 hold the thresholds of the step about to be decided (see make_thresholds)
    Rng rng[G];
    int32_t hf[G][NSL], ntk[G][NSL], l1[G][NSL];
    uint32_t sb[G], dif[G];
    int32_t flips_loc[G], nlin_cur[G];
#pragma unroll
    for (int g = 0; g < G; ++g) {
        const uint32_t trial = p.trial_offset + uint32_t((cl * G + g) * 4 + c4);
        rng_seed(rng[g], p.seed, trial, uint32_t(lane_l));
        uint32_t d[NSL]; draw_slots(rng[g], d);
        uint32_t s = ALLS;   // padding slots are +1
#pragma unroll
        for (int q = 0; q < NSL; ++q)
            if ((validm >> q) & 1) s = (s & ~(1u << q)) | ((d[q] >> 31) << q);
        sb[g] = s; dif[g] = s ^ ALLS;   // sprev = +1 everywhere at step 0
        flips_loc[g] = 0; nlin_cur[g] = 0;
#pragma unroll
        for (int q = 0; q < NSL; ++q) { hf[g][q] = b0[q]; ntk[g][q] = 0; l1[g][q] = 0; }
    }
    __syncthreads();
    cluster_sync_all();   // barriers of every CTA initialised before any st.async

    // remote addresses (fixed per thread). Lane L: destination CTA dst = L >> 2, role = L & 3.
    //   R = 2: lane sends the (F+, F) words of chain role (8 B at word wg*8 + 2 role).
    //   R = 1: role = (h, pp) sends the masks of chains {2pp, 2pp + 1} of B word W = 2 wg + h (16 B at W*8 + 4 pp).
    //   Role 0 also sends the warp's packed n_lin word (LIN builds).
    const int dst = L >> 2, role = L & 3;
    const bool sender = dst < CS;
    const bool msk_sender = sender;
    const bool cnt_sender = LIN && sender && role == 0;
    const int moff = (R == 2) ? wg * 8 + 2 * role : (2 * wg + (role >> 1)) * 8 + 4 * (role & 1);
    uint32_t r_msk[G][2], r_cnt[G][2], r_bar[G][2];
#pragma unroll
    for (int g = 0; g < G; ++g)
#pragma unroll
        for (int par = 0; par < 2; ++par) {
            const uint32_t base = smem_u32(RB + (g * 2 + par) * MW);
            const uint32_t d = sender ? uint32_t(dst) : 0u;
            r_msk[g][par] = mapa(base + 4u * uint32_t(moff), d);
            r_cnt[g][par] = mapa(base + 4u * uint32_t(MASKW + wg), d);
            r_bar[g][par] = mapa(bar0 + 8u * uint32_t(g * 2 + par), d);
        }

    // ---- send: assemble this warp's mask words (byte = RNG lane offset within the word's 4 lanes, bit = round) ----
    auto send = [&](int g, int m, uint32_t fbits, uint32_t snew, uint32_t lin) {
        const int par = m & 1;
        uint32_t cw = 0;
        if (LIN) cw = __reduce_add_sync(FULL, lin << (8 * c4));   // byte c = n_lin of chain c over this warp's slots (<= 64)
        uint32_t x0, y0, x1, y1;   // (F+, F) of the chain(s) this lane sends
        if constexpr (R == 1) {    // gather the 8 raw lanes (4 RNG lanes x 2 chains) of this lane's word and chain pair
            const uint32_t av = fbits | ((fbits & snew) << 8);
            const int h = role >> 1, pp = role & 1;
            uint32_t v[8];
#pragma unroll
            for (int i = 0; i < 8; ++i) v[i] = __shfl_sync(FULL, av, 4 * (4 * h + (i & 3)) + 2 * pp + (i >> 2));
            const uint32_t x01 = __byte_perm(v[0], v[1], 0x5140u), x23 = __byte_perm(v[2], v[3], 0x5140u);
            const uint32_t x45 = __byte_perm(v[4], v[5], 0x5140u), x67 = __byte_perm(v[6], v[7], 0x5140u);
            y0 = __byte_perm(x01, x23, 0x5410u); x0 = __byte_perm(x01, x23, 0x7632u);
            y1 = __byte_perm(x45, x67, 0x5410u); x1 = __byte_perm(x45, x67, 0x7632u);
        } else {                   // R = 2: each lane holds 4 rounds (2q + hsel) of RNG lane lam; disjoint bits -> add-reduce
            const uint32_t sh = uint32_t(8 * lam + hsel);
            uint32_t a = spread4(fbits) << sh, b = spread4(fbits & snew) << sh;
#pragma unroll
            for (int o = 4; o <= 16; o <<= 1) { a += __shfl_xor_sync(FULL, a, o); b += __shfl_xor_sync(FULL, b, o); }
            x0 = b; y0 = a; x1 = y1 = 0u;   // this lane's chain (role == c4)
        }
        const uint32_t rbr = par ? r_bar[g][1] : r_bar[g][0];
        if (R == 1) { if (msk_sender) st_async16(par ? r_msk[g][1] : r_msk[g][0], rbr, make_uint4(x0, y0, x1, y1)); }
        else if (msk_sender) st_async8(par ? r_msk[g][1] : r_msk[g][0], rbr, x0, y0);
        if (cnt_sender) st_async4(par ? r_cnt[g][1] : r_cnt[g][0], rbr, cw);
    };

    // thresholds of step t: ntk = q - thr (+/- corr folded for non-Onsager modes, using dif = flips of step t-1), l1 for the
    // n_lin band. The RNG stream does not depend on the chain state.
    auto make_thresholds = [&](int g, int t) {
        const int32_t fourT = T_fourT[t], twoT = fourT >> 1, qc = T_q[t], corr = T_kconst[t];
        const int32_t qs = qc - corr, qd = qc + corr;
        const int64_t c64 = -int64_t(32768) * fourT;
        const uint32_t df = dif[g];
        uint32_t d[NSL]; draw_slots(rng[g], d);
#pragma unroll
        for (int q = 0; q < NSL; ++q) {
            const int32_t thr = int32_t((int64_t(uint64_t(d[q] >> 16) * uint64_t(uint32_t(fourT))) + c64) >> 16);
            ntk[g][q] = ONS ? qc - thr : (((df >> q) & 1u) ? qd : qs) - thr;
            if (LIN) l1[g][q] = thr + twoT - 1;
        }
    };

    // ---- receive: wait for message m of group g, BMMAs, field update (factor 1 for the initial field, 2 afterwards) ----
    auto receive = [&](int g, int m) {
        const int par = m & 1;
        const uint32_t bar = bar0 + 8 * (g * 2 + par);
        const uint32_t parity = uint32_t(m >> 1) & 1u;
        while (!mbar_test(bar, parity)) {}
        PROF(5);
        if (tid == 0) mbar_arm(bar, BYTES);   // for message m + 2
        const uint32_t* rb = RB + (g * 2 + par) * MW;
        if (LIN) {   // n_lin: each lane unpacks NW/32 count words into 16-bit fields (chains 0/2 and 1/3), then two warp sums
            uint32_t s02 = 0, s13 = 0;
#pragma unroll
            for (int i = 0; i < NW / 32; ++i) { const uint32_t wv = rb[MASKW + L + 32 * i]; s02 += wv & 0x00FF00FFu; s13 += (wv >> 8) & 0x00FF00FFu; }
            s02 = __reduce_add_sync(FULL, s02); s13 = __reduce_add_sync(FULL, s13);
            const uint32_t sel = (c4 & 1) ? s13 : s02;
            nlin_cur[g] = int32_t((c4 & 2) ? (sel >> 16) : (sel & 0xFFFFu));
        }
        uint32_t b[8][2];
#pragma unroll
        for (int kb = 0; kb < 8; ++kb)
#pragma unroll
            for (int hi = 0; hi < 2; ++hi) b[kb][hi] = rb[(8 * kb + 4 * hi + (L & 3)) * 8 + (L >> 2)];
        int32_t acc[MAGS][2][JB][KS][4];
#pragma unroll
        for (int bb = 0; bb < MAGS; ++bb)
#pragma unroll
            for (int kd = 0; kd < 2; ++kd)
#pragma unroll
                for (int j = 0; j < JB; ++j)
#pragma unroll
                    for (int ks = 0; ks < KS; ++ks) acc[bb][kd][j][ks][0] = acc[bb][kd][j][ks][1] = acc[bb][kd][j][ks][2] = acc[bb][kd][j][ks][3] = 0;
#pragma unroll
        for (int kb = 0; kb < 8; ++kb)
#pragma unroll
            for (int j = 0; j < JB; ++j) {
                const uint4 sg = plane(0, j, kb);
#pragma unroll
                for (int bb = 0; bb < MAGS; ++bb) {
                    const uint4 mg = plane(1 + bb, j, kb);
                    // K = 2: |J| <= 1, so sg implies m_0 and m_0 & sg == sg exactly (no AND needed)
                    const uint4 ng = (MAGS == 1) ? sg : make_uint4(mg.x & sg.x, mg.y & sg.y, mg.z & sg.z, mg.w & sg.w);
                    bmma(acc[bb][0][j][kb % KS], mg, b[kb][0], b[kb][1]);
                    bmma(acc[bb][1][j][kb % KS], ng, b[kb][0], b[kb][1]);
                }
            }
        if (C::RNG_IN_RECV && m < S) make_thresholds(g, m);   // independent of the BMMAs: fills the tensor-pipe shadow
        const int32_t f = m == 0 ? 1 : 2;
#pragma unroll
        for (int q = 0; q < NSL; ++q) {
            const int j = q >> 1, hi = q & 1;
            int32_t dh = 0;
#pragma unroll
            for (int bb = 0; bb < MAGS; ++bb) {
                int32_t a1 = 0, b1 = 0, a2 = 0, b2 = 0;
#pragma unroll
                for (int ks = 0; ks < KS; ++ks) {
                    a1 += acc[bb][0][j][ks][2 * hi]; b1 += acc[bb][0][j][ks][2 * hi + 1];
                    a2 += acc[bb][1][j][ks][2 * hi]; b2 += acc[bb][1][j][ks][2 * hi + 1];
                }
                dh += (2 * a1 - b1 - 4 * a2 + 2 * b2) << bb;
            }
            hf[g][q] += f * dh;
        }
    };

#pragma unroll
    for (int g = 0; g < G; ++g) {
        if (!C::RNG_IN_RECV) make_thresholds(g, 0);
        send(g, 0, validm, sb[g], 0u);   // message 0: F := valid slots, sigma' := initial spins
    }

    const int chain_base = cl * G * 4;
    const bool writer = (wg == 0) && (L < 4);   // lane L = chain c4, RNG lane 0, round 0
    const int32_t hcl = WIDE ? p.hclamp : 0;
    for (int m = 0; m <= S; ++m) {
#pragma unroll
        for (int g = 0; g < G; ++g) {
            PROF(0);
            receive(g, m);
            PROF(1);
            const int chain = chain_base + g * 4 + c4;
            if (LIN && m >= 1 && p.do_trace && writer && chain < p.nchains) p.out_trace[size_t(chain) * S + (m - 1)] = nlin_cur[g];
            if (m == S) continue;
            const int t = m;
            // z - thr = s*h*65536 + q -/+ corr - thr  (corr subtracted if s == sprev, added if the spin flipped last step)
            int32_t ci = 0;
            if (ONS) ci = int32_t(((int64_t(nlin_cur[g]) * T_kcorr[t]) >> 8) + int64_t(T_kconst[t]));
            const uint32_t lw = uint32_t(2 * (T_fourT[t] >> 1) - 1);
            uint32_t fbits = 0, lin = 0;
            const uint32_t s = sb[g], df = dif[g];
#pragma unroll
            for (int q = 0; q < NSL; ++q) {
                int32_t hs = ((s >> q) & 1u) ? hf[g][q] : -hf[g][q];
                if (WIDE) hs = max(-hcl, min(hs, hcl));   // exact: decisions at |s h| >= hcl are sign-determined
                int32_t x = ntk[g][q];
                if (ONS) x += ((df >> q) & 1u) ? ci : -ci;
                const int32_t dv = hs * 65536 + x;   // z - thr (int32 range proven on the host)
                fbits |= (dv < 0 ? 1u : 0u) << q;
                if (LIN && ((validm >> q) & 1u)) lin += (uint32_t(dv + l1[g][q]) < lw) ? 1u : 0u;
            }
            fbits &= validm;
            sb[g] = s ^ fbits; dif[g] = fbits; flips_loc[g] += __popc(fbits);
            PROF(2);
            send(g, m + 1, fbits, sb[g], lin);
            PROF(3);
            if (!C::RNG_IN_RECV && t + 1 < S) make_thresholds(g, t + 1);   // overlaps the message flight
            PROF(4);
        }
    }
#ifdef SCA_PROFILE
    if (tid == 0 && p.prof) for (int i = 0; i < 8; ++i) p.prof[blockIdx.x * 8 + i] = pr[i];
#endif
#pragma unroll
    for (int g = 0; g < G; ++g) {
        const int chain = chain_base + g * 4 + c4;
        if (chain < p.nchains) {
#pragma unroll
            for (int q = 0; q < NSL; ++q) {
                const int slot = (R == 1 ? q : 2 * q + hsel) * 256 + lane_l;
                p.out_spins[size_t(chain) * NPAD + slot] = uint8_t((sb[g] >> q) & 1u);
                p.out_h[size_t(chain) * NPAD + slot] = hf[g][q];
            }
            p.out_flips[(size_t(chain) * NW + wg) * 8 + (L >> 2)] = flips_loc[g];
        }
    }
    cluster_sync_all();
}

}  // namespace kern

// ------------------------------------------------------------------------------------------------------------------------------
namespace host {

struct Device {
    int id = 0;
    cudaDeviceProp prop{};
    std::string uuid;
};

std::string uuid_str(const cudaUUID_t& u) {
    char b[64]; const unsigned char* x = reinterpret_cast<const unsigned char*>(u.bytes);
    std::snprintf(b, sizeof b, "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x", x[0], x[1], x[2], x[3], x[4], x[5],
                  x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    return b;
}

Device open_device() {
    int n = 0; CK(cudaGetDeviceCount(&n));
    if (n != 1) throw std::runtime_error("expected exactly one visible GPU (CUDA_VISIBLE_DEVICES=GPU-" + std::string(kExpectedUuid) + ")");
    Device d; CK(cudaGetDeviceProperties(&d.prop, 0)); d.uuid = uuid_str(d.prop.uuid);
    if (d.uuid != kExpectedUuid) throw std::runtime_error("wrong GPU " + d.uuid);
    CK(cudaSetDevice(0));
    return d;
}

// A fragments (same for every cluster): index ((((wg*K + pl)*JB + j)*8 + kb)*32 + L)*4 + r.
// Row rho of block j -> field slot y = round*256 + LPW*wg + rho % LPW with round = 2R*j + rho / LPW.
// k index kap of block kb -> word W = 8kb + kap/32 (= RNG lanes 4W..4W+3), x = 256*(kap%8) + 4W + (kap%32)/8.
// Bit = plane(J[x][y]) as run_trial reads it (h[y] += dense[x*NP + y] * s[x]); plane 0 = sign, plane 1+b = bit b of |J|.
std::vector<uint32_t> build_A(const jmat::IntJ& J, int K, int R) {
    const int NW = 32 * R, JB = 4 / R, LPW = 8 / R, M = (1 << (K - 1)) - 1;
    if (J.maxabs > M) throw std::runtime_error("|J| exceeds 2^(K-1)-1 for K=" + std::to_string(K));
    std::vector<uint32_t> A(size_t(NW) * K * JB * 8 * 32 * 4, 0);
    for (int wg = 0; wg < NW; ++wg)
        for (int pl = 0; pl < K; ++pl)
            for (int j = 0; j < JB; ++j)
                for (int kb = 0; kb < 8; ++kb)
                    for (int L = 0; L < 32; ++L)
                        for (int r = 0; r < 4; ++r) {
                            const int rho = (L >> 2) + ((r & 1) ? 8 : 0);
                            const int y = (2 * R * j + rho / LPW) * 256 + LPW * wg + rho % LPW;
                            const int kbase = ((r >> 1) ? 128 : 0) + 32 * (L & 3);
                            uint32_t word = 0;
                            for (int i = 0; i < 32; ++i) {
                                const int kap = kbase + i, W = 8 * kb + kap / 32, beta = kap % 32;
                                const int x = 256 * (beta % 8) + 4 * W + beta / 8;
                                const int v = (x < J.N && y < J.N) ? J.at(x, y) : 0;
                                const uint32_t bit = pl == 0 ? uint32_t(v < 0) : uint32_t((std::abs(v) >> (pl - 1)) & 1);
                                word |= bit << i;
                            }
                            A[((((size_t(wg) * K + pl) * JB + j) * 8 + kb) * 32 + L) * 4 + r] = word;
                        }
    return A;
}

// Integer-range plan. The kernel evaluates dv = s*h*65536 + x with x = q -/+ corr - thr, and dv + l1 = z + 2T - 1, in int32.
// |h_y| <= H = max_y sum_x |J_xy| + max_i |b_i| at all times, and |x|, |q -/+ corr + 2T - 1| <= xm = |q| + corr_max + 4T.
//  * plain kernel: needs H*65536 + xm < 2^31 at every step (as the coupling-only guard, with the bias added to H);
//  * WIDE kernel: the decision uses hc = clamp(s*h, -HC, HC) with 65536*HC > xm_max. For |s h| >= HC the flip test
//    (dv < 0) and the n_lin band test (-2T < z < 2T) are already decided by the sign of s*h, at hc exactly as at s*h, so the
//    outcome is unchanged; it needs 65536*HC + xm_max < 2^31 instead of a bound on H.
struct RangePlan { bool plain_ok = false; int32_t hclamp = 0; int64_t H = 0, xm_max = 0; };
RangePlan plan_ranges(const sca::Tables& tb, int64_t hmax, int64_t bmax, int N) {
    const int64_t lim = (int64_t(1) << 31) - 1;
    RangePlan rp; rp.H = hmax + bmax;
    for (size_t t = 0; t < tb.fourT.size(); ++t) {
        if (tb.fourT[t] <= 0) throw std::runtime_error("fourT must be positive");
        if (tb.kcorr[t] < 0) throw std::runtime_error("negative kcorr not supported");
        const int64_t corr_max = ((int64_t(N) * tb.kcorr[t]) >> 8) + std::llabs(int64_t(tb.kconst[t])) + 1;
        rp.xm_max = std::max<int64_t>(rp.xm_max, std::llabs(int64_t(tb.q[t])) + corr_max + int64_t(tb.fourT[t]));
    }
    rp.plain_ok = rp.H * 65536 + rp.xm_max < lim;
    const int64_t hc = rp.xm_max / 65536 + 1;   // smallest HC with 65536*HC > xm_max
    if (65536 * hc + rp.xm_max >= lim) throw std::runtime_error("table range exceeds int32 kernel arithmetic even with clamping");
    rp.hclamp = int32_t(hc);
    return rp;
}

std::string gpu_status() {
    std::string cmd = std::string("nvidia-smi --id=GPU-") + kExpectedUuid +
                      " --query-gpu=utilization.gpu,memory.used,power.draw,module.power.draw.average,clocks.sm --format=csv,noheader 2>/dev/null";
    FILE* f = popen(cmd.c_str(), "r"); if (!f) return "?";
    char buf[512] = {0}; size_t n = fread(buf, 1, sizeof buf - 1, f); pclose(f); buf[n] = 0;
    std::string s(buf); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

struct Launcher {
    int K = 2, R = 2, preg = 2, cs = 8, groups = 1;
    bool ons = false, lin = true;
    void* fn = nullptr;
    int threads = 0, nw = 0;
    size_t smem = 0;
    int max_clusters = 0;
    bool nonportable = false;
    bool wide = false;
    template <int KK, int RR, int PR, int CSS, int GG, bool WW> void pick_mode() {
        if (ons) fn = reinterpret_cast<void*>(kern::mb_kernel<KK, RR, PR, CSS, GG, true, true, WW>);
        else if (lin) fn = reinterpret_cast<void*>(kern::mb_kernel<KK, RR, PR, CSS, GG, false, true, WW>);
        else fn = reinterpret_cast<void*>(kern::mb_kernel<KK, RR, PR, CSS, GG, false, false, WW>);
    }
    template <int KK, int RR, int PR, int CSS, int GG> bool try_pick() {
        if (K != KK || R != RR || preg != PR || cs != CSS || groups != GG) return false;
        if (wide) pick_mode<KK, RR, PR, CSS, GG, true>(); else pick_mode<KK, RR, PR, CSS, GG, false>();
        using Cc = kern::Cfg<KK, RR, PR, CSS, GG>;
        smem = Cc::TAB_OFF; threads = Cc::THREADS; nw = Cc::NW;
        return true;
    }
    void setup(int S) {
        // Supported variants (see Cfg): K = 2: R2/P2/CS8 (J in registers, 8 CTAs x 8 warps), G = 1..3; R1/P1/CS8 (one plane in
        // shared memory, v2 geometry 8 CTAs x 4 warps), G = 1. K = 4: R2/P2/CS8 (two planes in shared memory), G = 1.
        const bool ok = try_pick<2, 2, 2, 8, 1>() || try_pick<2, 2, 2, 8, 2>() || try_pick<2, 2, 2, 8, 3>() ||
                        try_pick<2, 2, 2, 8, 4>() || try_pick<2, 1, 1, 8, 1>() || try_pick<2, 1, 1, 8, 2>() || try_pick<4, 2, 2, 8, 1>();
        if (!ok) throw std::runtime_error("unsupported variant (K/R/PREG/CS/G)");
        smem += size_t((12 * S + 15) / 16) * 16 + size_t(8) * S;
        CK(cudaFuncSetAttribute(fn, cudaFuncAttributeMaxDynamicSharedMemorySize, int(smem)));
        if (cs > 8) { CK(cudaFuncSetAttribute(fn, cudaFuncAttributeNonPortableClusterSizeAllowed, 1)); nonportable = true; }
        cudaLaunchConfig_t cfg = {}; cudaLaunchAttribute at[1];
        at[0].id = cudaLaunchAttributeClusterDimension; at[0].val.clusterDim.x = cs; at[0].val.clusterDim.y = 1; at[0].val.clusterDim.z = 1;
        cfg.attrs = at; cfg.numAttrs = 1; cfg.blockDim = dim3(threads); cfg.dynamicSmemBytes = smem; cfg.gridDim = dim3(cs);
        CK(cudaOccupancyMaxActiveClusters(&max_clusters, fn, &cfg));
    }
    void launch(int clusters, const kern::Params& p, cudaStream_t st) {
        cudaLaunchConfig_t cfg = {}; cudaLaunchAttribute at[1];
        at[0].id = cudaLaunchAttributeClusterDimension; at[0].val.clusterDim.x = cs; at[0].val.clusterDim.y = 1; at[0].val.clusterDim.z = 1;
        cfg.attrs = at; cfg.numAttrs = 1; cfg.blockDim = dim3(threads); cfg.dynamicSmemBytes = smem; cfg.gridDim = dim3(cs * clusters);
        cfg.stream = st;
        void* args[] = {const_cast<kern::Params*>(&p)};
        CK(cudaLaunchKernelExC(&cfg, fn, args));
    }
};

}  // namespace host

int main(int argc, char** argv) {
    try {
        v80host::Schedule sc;
        std::string jkind = "graph", jpath = "data/K2000.bin", output, mode = "run";
        int chains = 1, batches = 1, trace = 0, clusters_arg = 0, warmup = 1; double burn_s = 0;
        int K = 2, R = 2, preg = 2, cs = 8, groups = 1; int64_t target = INT64_MIN, target_energy = INT64_MAX;
        std::string bias_path, wide_arg = "auto"; int check_fields = 0;
        uint64_t seed = 20261004; uint32_t offset = 0;
        for (int a = 1; a < argc; ++a) {
            std::string k = argv[a];
            if (v80host::parse_schedule_flag(sc, k, a, argc, argv)) continue;
            if (k == "--info") { mode = "info"; continue; }
            auto val = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++a]); };
            if (k == "--graph") { jkind = "graph"; jpath = val(); }
            else if (k == "--gset") { jkind = "gset"; jpath = val(); }
            else if (k == "--jint8") { jkind = "jint8"; jpath = val(); }
            else if (k == "--target") target = std::stoll(val());
            else if (k == "--target-energy") target_energy = std::stoll(val());
            else if (k == "--bias") bias_path = val();
            else if (k == "--wide") wide_arg = val();
            else if (k == "--check-fields") check_fields = std::stoi(val());
            else if (k == "--kbits") K = std::stoi(val());
            else if (k == "--r") R = std::stoi(val());
            else if (k == "--preg") preg = std::stoi(val());
            else if (k == "--output") output = val();
            else if (k == "--chains") chains = std::stoi(val());
            else if (k == "--clusters") clusters_arg = std::stoi(val());
            else if (k == "--batches") batches = std::stoi(val());
            else if (k == "--cs") cs = std::stoi(val());
            else if (k == "--groups") groups = std::stoi(val());
            else if (k == "--trace") trace = std::stoi(val());
            else if (k == "--warmup") warmup = std::stoi(val());
            else if (k == "--burn") burn_s = std::stod(val());
            else if (k == "--seed") seed = std::stoull(val());
            else if (k == "--offset") offset = uint32_t(std::stoul(val()));
            else throw std::runtime_error("unknown argument " + k);
        }
        host::Device dev = host::open_device();
        const int S = sc.S;
        if (S < 2 || S > 4096) throw std::runtime_error("steps out of range");
        const jmat::IntJ J = jmat::load(jkind, jpath);
        const std::vector<int32_t> bias = bias_path.empty() ? std::vector<int32_t>(2048, 0) : jmat::load_bias(bias_path, J.N);
        const bool has_bias = !bias_path.empty();
        int64_t bmax = 0, bsum_abs = 0; for (int i = 0; i < 2048; ++i) { bmax = std::max<int64_t>(bmax, std::abs(bias[i])); bsum_abs += std::abs(bias[i]); }
        // success: --target CUT (Max-Cut, b = 0; default 33000 for K2000) or --target-energy E (energy <= E)
        if (target == INT64_MIN && target_energy == INT64_MAX) {
            if (jkind == "graph" && !has_bias) target = 33000; else throw std::runtime_error("--target CUT or --target-energy E required for this instance");
        }
        if (target != INT64_MIN && bsum_abs != 0) throw std::runtime_error("--target (cut) is defined for b = 0 only; use --target-energy");
        sca::Tables tb = v80host::build_tables(sc);
        const host::RangePlan rp = host::plan_ranges(tb, J.hmax, bmax, J.N);
        bool wide = false;
        if (wide_arg == "auto") wide = !rp.plain_ok;
        else if (wide_arg == "1") wide = true;
        else if (wide_arg == "0") { if (!rp.plain_ok) throw std::runtime_error("instance needs the WIDE kernel (|J| row sum + |b| too large for int32 decisions)"); }
        else throw std::runtime_error("--wide auto|0|1");
        bool ons = false; for (int t = 0; t < S; ++t) ons |= tb.kcorr[t] != 0;
        host::Launcher ln; ln.K = K; ln.R = R; ln.preg = preg; ln.cs = cs; ln.groups = groups; ln.ons = ons; ln.lin = ons || trace; ln.wide = wide; ln.setup(S);
        const int per_cluster = 4 * groups;
        if (mode == "info") {
            std::cout << "{\"device\":\"" << dev.prop.name << "\",\"uuid\":\"" << dev.uuid << "\",\"sms\":" << dev.prop.multiProcessorCount
                      << ",\"kbits\":" << K << ",\"r\":" << R << ",\"preg\":" << preg << ",\"onsager_kernel\":" << (ons ? "true" : "false")
                      << ",\"nlin_kernel\":" << (ln.lin ? "true" : "false") << ",\"cs\":" << cs << ",\"groups\":" << groups << ",\"threads\":" << ln.threads
                      << ",\"smem\":" << ln.smem << ",\"max_active_clusters\":" << ln.max_clusters << ",\"chains_per_wave\":" << ln.max_clusters * per_cluster
                      << ",\"N\":" << J.N << ",\"maxabs\":" << J.maxabs << ",\"hmax\":" << J.hmax << ",\"bmax\":" << bmax << ",\"H\":" << rp.H
                      << ",\"plain_ok\":" << (rp.plain_ok ? "true" : "false") << ",\"wide\":" << (wide ? "true" : "false") << ",\"hclamp\":" << rp.hclamp
                      << ",\"sumw\":" << J.sumw << "}\n";
            return 0;
        }
        if (output.empty()) throw std::runtime_error("--output PREFIX required");
        if (burn_s > 0 && std::ifstream(output + ".burn.json").good()) throw std::runtime_error("refusing to overwrite " + output);
        if (std::ifstream(output + ".trials.jsonl").good()) throw std::runtime_error("refusing to overwrite " + output);
        const int clusters = clusters_arg > 0 ? clusters_arg : (chains + per_cluster - 1) / per_cluster;
        if (chains < 1 || chains > clusters * per_cluster) throw std::runtime_error("chains must be in 1..clusters*4*groups");

        const std::vector<uint32_t> A = host::build_A(J, K, R);
        const int NW = ln.nw;
        uint4* dA; int32_t *dFourT, *dQ, *dKconst, *dH, *dFl, *dTrace = nullptr; int64_t* dKcorr; uint8_t* dSpins;
        CK(cudaMalloc(&dA, A.size() * 4)); CK(cudaMemcpy(dA, A.data(), A.size() * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dFourT, S * 4)); CK(cudaMemcpy(dFourT, tb.fourT.data(), S * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dQ, S * 4)); CK(cudaMemcpy(dQ, tb.q.data(), S * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dKconst, S * 4)); CK(cudaMemcpy(dKconst, tb.kconst.data(), S * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dKcorr, S * 8)); CK(cudaMemcpy(dKcorr, tb.kcorr.data(), S * 8, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dSpins, size_t(chains) * 2048)); CK(cudaMalloc(&dH, size_t(chains) * 2048 * 4));
        CK(cudaMalloc(&dFl, size_t(chains) * NW * 8 * 4));
        if (trace) CK(cudaMalloc(&dTrace, size_t(chains) * S * 4));
        int32_t* dBias = nullptr;
        if (has_bias) { CK(cudaMalloc(&dBias, 2048 * 4)); CK(cudaMemcpy(dBias, bias.data(), 2048 * 4, cudaMemcpyHostToDevice)); }
        kern::Params p{};
        p.A = dA; p.fourT = dFourT; p.q = dQ; p.kcorr = dKcorr; p.kconst = dKconst; p.S = S; p.N = J.N; p.seed = seed; p.nchains = chains;
        p.do_trace = trace; p.out_spins = dSpins; p.out_h = dH; p.out_flips = dFl; p.out_trace = dTrace;
        p.bias = dBias; p.hclamp = rp.hclamp;
#ifdef SCA_PROFILE
        long long* dProf; CK(cudaMalloc(&dProf, size_t(clusters) * cs * 8 * 8)); CK(cudaMemset(dProf, 0, size_t(clusters) * cs * 8 * 8)); p.prof = dProf;
#endif
        cudaStream_t st; CK(cudaStreamCreateWithFlags(&st, cudaStreamNonBlocking));
        cudaEvent_t e0, e1; CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));
        // warm-up launches on trial ids far outside every cohort (results discarded)
        for (int i = 0; i < warmup; ++i) { p.trial_offset = 0xF0000000u + uint32_t(i) * 65536u; ln.launch(clusters, p, st); }
        CK(cudaStreamSynchronize(st));

        if (burn_s > 0) {   // sustained back-to-back batches for power sampling; outcomes are not read (not samples)
            std::ofstream bj(output + ".burn.json");
            const auto b0 = std::chrono::steady_clock::now(); const std::string st0 = host::gpu_status();
            long long nb = 0; double dsum = 0;
            CK(cudaEventRecord(e0, st));
            while (std::chrono::duration<double>(std::chrono::steady_clock::now() - b0).count() < burn_s) {
                for (int i = 0; i < 64; ++i) { p.trial_offset = offset + uint32_t((nb % 4096) * chains); ln.launch(clusters, p, st); ++nb; }
                CK(cudaStreamSynchronize(st));
            }
            CK(cudaEventRecord(e1, st)); CK(cudaEventSynchronize(e1)); float ms = 0; CK(cudaEventElapsedTime(&ms, e0, e1)); dsum = ms;
            const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - b0).count();
            bj << std::setprecision(10) << "{\"burn_seconds_wall\":" << wall << ",\"batches\":" << nb << ",\"chains\":" << chains << ",\"kbits\":" << K
               << ",\"instance\":\"" << jpath << "\",\"bias\":\"" << bias_path << "\",\"wide\":" << (wide ? "true" : "false")
               << ",\"r\":" << R << ",\"preg\":" << preg << ",\"cs\":" << cs << ",\"groups\":" << groups << ",\"steps\":" << S
               << ",\"device_ms_total\":" << dsum << ",\"device_ms_per_batch\":" << dsum / nb
               << ",\"gpu_before\":\"" << st0 << "\",\"gpu_after\":\"" << host::gpu_status() << "\"}\n";
            std::cout << "burn: " << nb << " batches in " << wall << " s, " << dsum / nb << " ms/batch\n";
            return 0;
        }
        std::ofstream trials(output + ".trials.jsonl"), bat(output + ".batches.jsonl"), spins(output + ".spins.bin", std::ios::binary);
        std::ofstream traces;
        if (trace) traces.open(output + ".trace.bin", std::ios::binary);
        if (!trials || !bat || !spins) throw std::runtime_error("cannot create output files");
        const std::string status_before = host::gpu_status();
        std::vector<uint8_t> hs(size_t(chains) * 2048); std::vector<int32_t> hh(size_t(chains) * 2048), hfl(size_t(chains) * NW * 8);
        std::vector<int32_t> htr(trace ? size_t(chains) * S : 0);
        int64_t succ = 0, batch_succ = 0, mism = 0; double dev_ms_sum = 0, host_ms_sum = 0, dev_ms_min = 1e30, dev_ms_max = 0;
        for (int b = 0; b < batches; ++b) {
            p.trial_offset = offset + uint32_t(b) * uint32_t(chains);
            const auto h0 = std::chrono::steady_clock::now();
            CK(cudaEventRecord(e0, st));
            ln.launch(clusters, p, st);
            CK(cudaEventRecord(e1, st));
            CK(cudaEventSynchronize(e1));
            const auto h1 = std::chrono::steady_clock::now();
            float ms = 0; CK(cudaEventElapsedTime(&ms, e0, e1));
            const double hms = std::chrono::duration<double, std::milli>(h1 - h0).count();
            CK(cudaGetLastError());
            CK(cudaMemcpy(hs.data(), dSpins, hs.size(), cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(hh.data(), dH, hh.size() * 4, cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(hfl.data(), dFl, hfl.size() * 4, cudaMemcpyDeviceToHost));
            if (trace) CK(cudaMemcpy(htr.data(), dTrace, htr.size() * 4, cudaMemcpyDeviceToHost));
            int bs = 0;
            std::vector<std::array<uint64_t, 32>> bitsv(chains);
            std::vector<int64_t> cdev(chains), chost(chains), flips(chains), shv(chains), bsv(chains), edev(chains), ehost(chains), fmis(chains, 0);
            {
                std::atomic<int> nx{0};
                auto work = [&]() {
                    for (int c; (c = nx++) < chains;) {
                        uint64_t* bits = bitsv[c].data(); std::fill(bits, bits + 32, 0ull); int64_t sh = 0, sbs = 0;
                        for (int i = 0; i < J.N; ++i) {
                            const bool up = hs[size_t(c) * 2048 + i] != 0;
                            if (up) bits[i >> 6] |= uint64_t(1) << (i & 63);
                            const int64_t h = hh[size_t(c) * 2048 + i];
                            sh += up ? h : -h;
                            sbs += up ? bias[i] : -bias[i];
                        }
                        shv[c] = sh; bsv[c] = sbs;
                        edev[c] = -(sh + sbs) / 2;                      // device energy E = -(sum s h + sum b s) / 2
                        ehost[c] = jmat::energy(J, bias.data(), bits);  // independent: -sum_{i<j} J s s - sum b s
                        cdev[c] = (J.sumw + sh / 2) / 2;                // Max-Cut value (b = 0), as run_trial
                        chost[c] = bsum_abs == 0 ? jmat::scalar_cut(J, bits) : 0;
                        if (check_fields) fmis[c] = jmat::field_mismatches(J, bias.data(), bits, &hh[size_t(c) * 2048]);
                        int64_t f = 0; for (int i = 0; i < NW * 8; ++i) f += hfl[size_t(c) * NW * 8 + i];
                        flips[c] = f;
                    }
                };
                const int nth = std::min(chains, 48);
                std::vector<std::thread> th; for (int i = 0; i < nth; ++i) th.emplace_back(work); for (auto& x : th) x.join();
            }
            for (int c = 0; c < chains; ++c) {
                const int64_t cut_dev = cdev[c], cut_host = chost[c];
                if (edev[c] != ehost[c] || (bsum_abs == 0 && cut_dev != cut_host) || fmis[c] != 0) ++mism;
                const bool ok = target != INT64_MIN ? cut_host >= target : ehost[c] <= target_energy;
                succ += ok; bs += ok;
                const uint32_t tid = p.trial_offset + uint32_t(c);
                trials << "{\"batch\":" << b << ",\"chain\":" << c << ",\"trial_id\":" << tid;
                if (bsum_abs == 0) trials << ",\"cut\":" << cut_host << ",\"device_cut\":" << cut_dev;
                trials << ",\"sum_sh\":" << shv[c] << ",\"sum_bs\":" << bsv[c] << ",\"energy\":" << ehost[c] << ",\"device_energy\":" << edev[c];
                if (check_fields) trials << ",\"field_mismatches\":" << fmis[c];
                trials << ",\"flips\":" << flips[c] << ",\"success\":" << (ok ? "true" : "false") << "}\n";
                spins.write(reinterpret_cast<const char*>(bitsv[c].data()), 32 * 8);
                if (trace) traces.write(reinterpret_cast<const char*>(&htr[size_t(c) * S]), size_t(S) * 4);
            }
            batch_succ += bs > 0;
            dev_ms_sum += ms; host_ms_sum += hms; dev_ms_min = std::min<double>(dev_ms_min, ms); dev_ms_max = std::max<double>(dev_ms_max, ms);
            bat << std::setprecision(9) << "{\"batch\":" << b << ",\"chains\":" << chains << ",\"trial_offset\":" << p.trial_offset
                << ",\"device_ms\":" << ms << ",\"host_ms\":" << hms << ",\"successes\":" << bs << "}\n";
        }
#ifdef SCA_PROFILE
        { std::vector<long long> pr(size_t(clusters) * cs * 8); CK(cudaMemcpy(pr.data(), dProf, pr.size() * 8, cudaMemcpyDeviceToHost));
          const char* nm[8] = {"loop", "recv_compute", "decide", "send", "thresholds", "wait", "-", "-"};
          for (int i = 0; i < 6; ++i) { double s = 0; for (int b = 0; b < clusters * cs; ++b) s += pr[size_t(b) * 8 + i];
            std::cout << "prof " << nm[i] << " cycles/step/group = " << s / (clusters * cs) / (double(S) * groups) << "\n"; } }
#endif
        const std::string status_after = host::gpu_status();
        const double tb_mean = dev_ms_sum / batches, pt = double(succ) / (double(chains) * batches), Pb = double(batch_succ) / batches;
        auto tts = [&](double P, double t) { return P >= 1.0 ? t : (P <= 0 ? -1.0 : t * std::log(0.01) / std::log(1.0 - P)); };
        auto wilson = [](double k, double n, int side) {   // 95% Wilson score interval
            const double z = 1.959963984540054, ph = k / n, den = 1 + z * z / n, cen = (ph + z * z / (2 * n)) / den;
            const double half = z * std::sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den;
            return side < 0 ? std::max(0.0, cen - half) : std::min(1.0, cen + half);
        };
        std::ofstream sum(output + ".summary.json");
        sum << std::setprecision(10) << "{\"device\":\"" << dev.prop.name << "\",\"uuid\":\"" << dev.uuid << "\",\"kbits\":" << K << ",\"r\":" << R
            << ",\"preg\":" << preg << ",\"cs\":" << cs << ",\"groups\":" << groups
            << ",\"clusters\":" << clusters << ",\"max_active_clusters\":" << ln.max_clusters << ",\"threads\":" << ln.threads << ",\"smem\":" << ln.smem
            << ",\"instance\":\"" << jpath << "\",\"instance_kind\":\"" << jkind << "\",\"N\":" << J.N << ",\"sumw\":" << J.sumw
            << ",\"bias\":\"" << bias_path << "\",\"bmax\":" << bmax << ",\"hmax\":" << J.hmax << ",\"wide\":" << (wide ? "true" : "false") << ",\"hclamp\":" << rp.hclamp
            << ",\"plain_ok\":" << (rp.plain_ok ? "true" : "false") << ",\"check_fields\":" << check_fields
            << ",\"target\":" << (target == INT64_MIN ? std::string("null") : std::to_string(target))
            << ",\"target_energy\":" << (target_energy == INT64_MAX ? std::string("null") : std::to_string(target_energy))
            << ",\"chains\":" << chains << ",\"batches\":" << batches << ",\"seed\":" << seed << ",\"trial_offset\":" << offset
            << ",\"onsager_kernel\":" << (ons ? "true" : "false") << ",\"nlin_kernel\":" << (ln.lin ? "true" : "false")
            << ",\"steps\":" << S << ",\"t0\":" << sc.t0 << ",\"t1\":" << sc.t1 << ",\"q\":" << sc.q << ",\"lambda\":" << sc.lam
            << ",\"ramp\":" << (sc.ramp ? "true" : "false") << ",\"tec_jv\":" << sc.jv << ",\"tecT_kappa\":" << sc.kappa << ",\"trace\":" << trace
            << ",\"successes\":" << succ << ",\"trials\":" << int64_t(chains) * batches << ",\"p\":" << pt << ",\"batches_with_success\":" << batch_succ
            << ",\"P_batch\":" << Pb << ",\"p_wilson95\":[" << wilson(double(succ), double(chains) * batches, -1) << "," << wilson(double(succ), double(chains) * batches, 1) << "]"
            << ",\"P_batch_wilson95\":[" << wilson(double(batch_succ), batches, -1) << "," << wilson(double(batch_succ), batches, 1) << "]"
            << ",\"device_ms_mean\":" << tb_mean << ",\"device_ms_min\":" << dev_ms_min << ",\"device_ms_max\":" << dev_ms_max
            << ",\"host_ms_mean\":" << host_ms_sum / batches << ",\"tts99_primary_ms\":" << tts(Pb, tb_mean)
            << ",\"cut_mismatches\":" << mism << ",\"mismatches\":" << mism
            << ",\"gpu_before\":\"" << status_before << "\",\"gpu_after\":\"" << status_after << "\"}\n";
        std::cout << "chains=" << chains << " batches=" << batches << " p=" << pt << " P_batch=" << Pb << " t_batch=" << tb_mean
                  << " ms  TTS_primary=" << tts(Pb, tb_mean) << " ms  wide=" << (wide ? 1 : 0) << " mismatch=" << mism << "\n";
        return mism ? 2 : 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
