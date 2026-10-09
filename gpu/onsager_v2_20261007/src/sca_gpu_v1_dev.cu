// Bit-exact GPU (GH200, sm_90) implementation of the V80 Onsager/TEC-T/plain SCA engine (reference: sca_ref.hpp, SCA_LANES=256).
//
// One thread-block cluster of CS CTAs runs 4*G independent chains. The cluster holds the whole coupling matrix J in registers as
// binary tensor-core (BMMA m16n8k256 .b1 AND+POPC) A-fragments: 32 warps per cluster, each owning 8 RNG lanes x 8 rounds = 64
// field slots y and the full 2048-bit x range (128 registers per thread). Lane L of a warp owns chain (L % 4) of each group and
// RNG lane 8*wg + L/4, i.e. the 8 slots k*256 + l (k = 0..7) that lane l of the V80 engine decides in rounds k = 0..7.
// Per step and group, every warp:
//   decisions (8 per thread, exactly as sca::run_trial, Q16.16 integer arithmetic) -> per-chain flip masks F and F&sigma'
//   -> CTA staging -> st.async (DSMEM, mbarrier complete_tx) to every CTA of the cluster
//   -> after the mbarrier wait, BMMA over the 8 k-blocks gives popc(J_y & Fsigma') and popc(J_y & F) for its 64 y and 4 chains
//   -> dh_y = 2*(4 cA - 2 cB - 2|Fsigma'| + |F| + [y in F](2 sigma'_y - 1))   (diagonal and padding handled exactly)
// The initial field uses the same path with F := valid slots, sigma' := initial spins and factor 1.
// Random numbers do not depend on the chain state; each step's 8 thresholds per thread are generated while messages are in flight.
#include "sca_ref.hpp"
#include "tables.hpp"
#include <cuda_runtime.h>
#include <algorithm>
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
#include <vector>

#if SCA_LANES != 256
#error "build with -DSCA_LANES=256"
#endif

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) throw std::runtime_error(std::string(#x) + ": " + cudaGetErrorString(e_)); } while (0)

static const char* kExpectedUuid = "8f88f1cd-2cd1-ab46-b984-2e97ebd18861";

namespace kern {

constexpr int N = 2000;

struct Params {
    const uint4* A;          // [32 wg][4 j][8 kb][32 lanes]
    const int32_t* fourT;    // [S]
    const int32_t* q;        // [S]
    const int64_t* kcorr;    // [S]
    const int32_t* kconst;   // [S]
    int S;
    uint64_t seed;
    uint32_t trial_offset;
    int nchains;             // chains reported (chain index < nchains); other chains of a group are computed but not written
    int do_trace;
    uint8_t* out_spins;      // [nchains][256]: bit k = (spin of slot k*256+l == +1)
    int32_t* out_sh;         // [nchains][256]: sum over the lane's valid slots of s*h
    int64_t* out_flips;      // [nchains]
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
__device__ __forceinline__ void st_async16(uint32_t raddr, uint32_t rbar, uint4 v) {
    asm volatile("st.async.shared::cluster.mbarrier::complete_tx::bytes.v4.b32 [%0], {%1,%2,%3,%4}, [%5];"
                 :: "r"(raddr), "r"(v.x), "r"(v.y), "r"(v.z), "r"(v.w), "r"(rbar) : "memory");
}
__device__ __forceinline__ void bmma(int32_t (&c)[4], const uint4& a, uint32_t b0, uint32_t b1) {
    asm("mma.sync.aligned.m16n8k256.row.col.s32.b1.b1.s32.and.popc {%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%0,%1,%2,%3};"
        : "+r"(c[0]), "+r"(c[1]), "+r"(c[2]), "+r"(c[3]) : "r"(a.x), "r"(a.y), "r"(a.z), "r"(a.w), "r"(b0), "r"(b1));
}

// xoshiro128**, seeded by splitmix64(seed, trial, lane) -- identical to sca::Rng.
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
// r = ((u16 - 32768) * fourT) >> 16 (int64, arithmetic shift); the result always fits in int32.
__device__ __forceinline__ int32_t threshold(uint32_t draw, int32_t fourT) {
    const int32_t a = int32_t(draw >> 16) - 32768;
    return int32_t((int64_t(a) * int64_t(fourT)) >> 16);
}

template <int CS, int G>
struct Cfg {
    static constexpr int NWG = 32 / CS;                 // warps per CTA
    static constexpr int THREADS = 32 * NWG;
    static constexpr int MSGW = NWG * 16 + 8;           // words per CTA message: masks [NWG][8 n][2 h'] + counts [4 chains] u64
    static constexpr int NCH = MSGW / 4;                // 16-byte chunks per message
    static constexpr uint32_t BYTES = uint32_t(CS * MSGW * 4);   // bytes expected per message (all CS senders)
    static constexpr int RB_OFF = 0;                                    // [G][2][CS][MSGW] u32
    static constexpr int ST_OFF = RB_OFF + G * 2 * CS * MSGW * 4;       // [G][NWG*16] u32
    static constexpr int CNT_OFF = ST_OFF + G * NWG * 16 * 4;           // [G][2][4] u64
    static constexpr int BAR_OFF = CNT_OFF + G * 2 * 4 * 8;             // [G][2] u64
    static constexpr int TAB_OFF = BAR_OFF + G * 2 * 8;                 // fourT[S], q[S], kconst[S] (i32), kcorr[S] (i64, 16-aligned)
    static_assert(MSGW % 4 == 0 && TAB_OFF % 16 == 0, "alignment");
};

// ONS: the correction depends on the previous step's n_lin (lambda > 0). Otherwise corr = kconst[t] is known in advance and is
// folded into the thresholds. LIN: count n_lin (needed for ONS and for the trace).
template <int CS, int G, bool ONS, bool LIN>
__global__ void __launch_bounds__(32 * (32 / CS), 1) sca_kernel(const Params p) {
    static_assert(LIN || !ONS, "Onsager mode needs n_lin");
    using C = Cfg<CS, G>;
    constexpr int NWG = C::NWG, MSGW = C::MSGW, NCH = C::NCH;
    extern __shared__ __align__(16) uint8_t smem[];
    uint32_t* RB = reinterpret_cast<uint32_t*>(smem + C::RB_OFF);
    uint32_t* ST = reinterpret_cast<uint32_t*>(smem + C::ST_OFF);
    unsigned long long* CNT = reinterpret_cast<unsigned long long*>(smem + C::CNT_OFF);
    const int S = p.S;
    int32_t* T_fourT = reinterpret_cast<int32_t*>(smem + C::TAB_OFF);
    int32_t* T_q = T_fourT + S;
    int32_t* T_kconst = T_q + S;
    int64_t* T_kcorr = reinterpret_cast<int64_t*>(smem + C::TAB_OFF + ((12 * S + 15) / 16) * 16);

    const int tid = threadIdx.x, L = tid & 31, w = tid >> 5;
    const uint32_t rank = cluster_rank();
    const int wg = int(rank) * NWG + w;                 // warp-global index 0..31
    const int c4 = L & 3;                               // chain within group
    const int lane_l = 8 * wg + (L >> 2);               // RNG lane 0..255
    const uint32_t validm = lane_l < (N - 7 * 256) ? 0xFFu : 0x7Fu;   // slot 7*256+l exists only for l < 208
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
    for (int i = tid; i < G * 8; i += blockDim.x) CNT[i] = 0ull;
    const uint32_t bar0 = smem_u32(smem + C::BAR_OFF);
    if (tid == 0) {
        for (int i = 0; i < 2 * G; ++i) mbar_init(bar0 + 8 * i, 1);
        asm volatile("fence.mbarrier_init.release.cluster;" ::: "memory");
        for (int i = 0; i < 2 * G; ++i) mbar_arm(bar0 + 8 * i, C::BYTES);   // messages 0 and 1 of every group
    }
    // A fragments: 4 y-blocks x 8 k-blocks x 4 registers, resident for the whole kernel.
    uint4 A[4][8];
#pragma unroll
    for (int j = 0; j < 4; ++j)
#pragma unroll
        for (int kb = 0; kb < 8; ++kb) A[j][kb] = p.A[((wg * 4 + j) * 8 + kb) * 32 + L];

    // per-group chain state. hh = h - s (shifted field: no diagonal correction needed, see header).
    Rng rng[G];
    int32_t hh[G][8], ntk[G][8], l1[G][8];
    uint32_t sb[G], dif[G];
    int32_t flips_tot[G], nlin_cur[G];
#pragma unroll
    for (int g = 0; g < G; ++g) {
        const uint32_t trial = p.trial_offset + uint32_t((cl * G + g) * 4 + c4);
        rng_seed(rng[g], p.seed, trial, uint32_t(lane_l));
        uint32_t s = 0xFFu;   // padding slots are +1
#pragma unroll
        for (int k = 0; k < 8; ++k) {
            const uint32_t r = rng_next(rng[g]);
            if ((validm >> k) & 1) s = (s & ~(1u << k)) | ((r >> 31) << k);
        }
        sb[g] = s; dif[g] = s ^ 0xFFu;   // sprev = +1 everywhere at step 0
        flips_tot[g] = 0; nlin_cur[g] = 0;
#pragma unroll
        for (int k = 0; k < 8; ++k) { hh[g][k] = 0; ntk[g][k] = 0; l1[g][k] = 0; }
    }
    __syncthreads();
    cluster_sync_all();   // barriers of every CTA initialised before any st.async

    auto boff = [&](int kb, int hi) -> int {   // B-fragment word for k-block kb (hi: second register)
        const int wgx = 4 * kb + 2 * hi + ((L >> 1) & 1);
        return (wgx / NWG) * MSGW + (wgx % NWG) * 16 + (L >> 2) * 2 + (L & 1);
    };

    // ---- send: pack this warp's per-chain masks (5 instructions), per-CTA counts (smem atomics), barrier, st.async ----
    auto send = [&](int g, int m, uint32_t fbits, uint32_t snew, uint32_t lin) {
        uint32_t* st = ST + g * (NWG * 16);
        const int par = m & 1;
        const uint32_t fs = fbits & snew;
        const uint32_t a = fbits | (fs << 16);
        const uint32_t b = __shfl_down_sync(0xffffffffu, a, 4);    // lane offset +1
        const uint32_t c = a | (b << 8);
        const uint32_t d = __shfl_down_sync(0xffffffffu, c, 8);    // lane offset +2
        if ((L & 12) == 0) {   // lane offsets 0 and 4: words for lanes [4h', 4h'+4) of this warp, chain c4
            const int hp = L >> 4;
            st[w * 16 + 4 * c4 + hp] = __byte_perm(c, d, 0x7632u);       // F & sigma'  (column 2c4)
            st[w * 16 + 4 * c4 + 2 + hp] = __byte_perm(c, d, 0x5410u);   // F           (column 2c4+1)
        }
        uint32_t pc = uint32_t(__popc(fbits)) | (uint32_t(__popc(fs)) << 10) | (lin << 20);
        pc += __shfl_xor_sync(0xffffffffu, pc, 4);
        pc += __shfl_xor_sync(0xffffffffu, pc, 8);
        pc += __shfl_xor_sync(0xffffffffu, pc, 16);
        unsigned long long* cnt = CNT + (g * 2 + par) * 4;
        if (L < 4) atomicAdd(cnt + L, (unsigned long long)(pc & 1023u) | ((unsigned long long)((pc >> 10) & 1023u) << 21) |
                                          ((unsigned long long)(pc >> 20) << 42));
        __syncthreads();
        if (tid == 0) { unsigned long long* o = CNT + (g * 2 + (par ^ 1)) * 4; o[0] = o[1] = o[2] = o[3] = 0ull; }
        const uint32_t my_rb = smem_u32(RB + ((g * 2 + par) * CS + int(rank)) * MSGW);
        const uint32_t my_bar = bar0 + 8 * (g * 2 + par);
#pragma unroll
        for (int it = 0; it < (CS * NCH + C::THREADS - 1) / C::THREADS; ++it) {
            const int tau = tid + it * C::THREADS;
            if (tau < CS * NCH) {
                const int dd = tau / NCH, qc = tau % NCH;
                const uint4 val = qc < NWG * 4 ? *reinterpret_cast<const uint4*>(st + 4 * qc)
                                               : *reinterpret_cast<const uint4*>(cnt + 2 * (qc - NWG * 4));
                st_async16(mapa(my_rb + 16u * qc, uint32_t(dd)), mapa(my_bar, uint32_t(dd)), val);
            }
        }
    };

    // ---- receive: wait for message m of group g, BMMA, field update (factor 1 for the initial field, 2 afterwards) ----
    auto receive = [&](int g, int m) {
        const int par = m & 1;
        const uint32_t bar = bar0 + 8 * (g * 2 + par);
        const uint32_t parity = uint32_t(m >> 1) & 1u;
        while (!mbar_test(bar, parity)) {}
        PROF(5);
        if (tid == 0) mbar_arm(bar, C::BYTES);   // for message m + 2
        const uint32_t* rb = RB + (g * 2 + par) * CS * MSGW;
        unsigned long long cs = 0;
#pragma unroll
        for (int s = 0; s < CS; ++s) cs += *reinterpret_cast<const unsigned long long*>(rb + s * MSGW + NWG * 16 + 2 * c4);
        const int32_t nF = int32_t(cs & 0x1FFFFFull), nFs = int32_t((cs >> 21) & 0x1FFFFFull);
        nlin_cur[g] = int32_t(cs >> 42);
        uint32_t b[8][2];
#pragma unroll
        for (int kb = 0; kb < 8; ++kb) { b[kb][0] = rb[boff(kb, 0)]; b[kb][1] = rb[boff(kb, 1)]; }
        int32_t acc[4][4];
#pragma unroll
        for (int j = 0; j < 4; ++j) {
            acc[j][0] = acc[j][1] = acc[j][2] = acc[j][3] = 0;
#pragma unroll
            for (int kb = 0; kb < 8; ++kb) bmma(acc[j], A[j][kb], b[kb][0], b[kb][1]);
        }
        const int32_t base = nF - 2 * nFs, f = m == 0 ? 1 : 2;
#pragma unroll
        for (int k = 0; k < 8; ++k) {
            const int j = k >> 1, hi = k & 1;
            hh[g][k] += f * (4 * acc[j][2 * hi] - 2 * acc[j][2 * hi + 1] + base);
        }
        if (m >= 1) flips_tot[g] += nF;
    };

    // thresholds for step t (random stream is state independent): ntk = q + 65536 - thr [- corr*e for non-Onsager], l1 for n_lin
    auto make_thresholds = [&](int g, int t) {
        const int32_t fourT = T_fourT[t], twoT = fourT >> 1, qc = T_q[t] + 65536;
        const int32_t corr = T_kconst[t];   // non-Onsager: corr = kconst[t] (kcorr == 0)
#pragma unroll
        for (int k = 0; k < 8; ++k) {
            const int32_t thr = threshold(rng_next(rng[g]), fourT);
            if (ONS) ntk[g][k] = qc - thr;
            else ntk[g][k] = qc - thr + (((dif[g] >> k) & 1u) ? corr : -corr);
            if (LIN) l1[g][k] = thr + twoT - 1;
        }
    };

#pragma unroll
    for (int g = 0; g < G; ++g) send(g, 0, validm, sb[g], 0u);   // message 0: F := valid slots, sigma' := initial spins
#pragma unroll
    for (int g = 0; g < G; ++g) make_thresholds(g, 0);

    const int chain_base = cl * G * 4;
    const bool writer = (wg == 0) && (L < 4);   // lane L = chain c4, RNG lane 0
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
            int32_t ci = 0;
            if (ONS) ci = int32_t(((int64_t(nlin_cur[g]) * T_kcorr[t]) >> 8) + int64_t(T_kconst[t]));
            const uint32_t lw = uint32_t(2 * (T_fourT[t] >> 1) - 1);
            uint32_t fbits = 0, lin = 0;
            const uint32_t s = sb[g], df = dif[g];
#pragma unroll
            for (int k = 0; k < 8; ++k) {
                const int32_t hs = ((s >> k) & 1u) ? hh[g][k] : -hh[g][k];
                int32_t x = ntk[g][k];
                if (ONS) x += ((df >> k) & 1u) ? ci : -ci;
                const int32_t dv = hs * 65536 + x;   // z - thr
                fbits |= (dv < 0 ? 1u : 0u) << k;
                if (LIN && (k < 7 || validm == 0xFFu)) lin += (uint32_t(dv + l1[g][k]) < lw) ? 1u : 0u;
            }
            fbits &= validm;
            sb[g] = s ^ fbits; dif[g] = fbits;
            PROF(2);
            send(g, m + 1, fbits, sb[g], lin);
            PROF(3);
            if (t + 1 < S) make_thresholds(g, t + 1);
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
            int32_t sh = 0;
#pragma unroll
            for (int k = 0; k < 8; ++k)
                if ((validm >> k) & 1u) sh += (((sb[g] >> k) & 1u) ? hh[g][k] : -hh[g][k]) + 1;   // s*h = s*hh + 1
            p.out_spins[size_t(chain) * 256 + lane_l] = uint8_t(sb[g]);
            p.out_sh[size_t(chain) * 256 + lane_l] = sh;
            if (writer) p.out_flips[chain] = flips_tot[g];
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

// A fragments for the whole matrix (same for every cluster): index ((wg*4 + j)*8 + kb)*32 + L, 4 words each.
std::vector<uint32_t> build_A(const sca::Graph& g) {
    std::vector<uint32_t> A(size_t(32) * 4 * 8 * 32 * 4, 0);
    for (int wg = 0; wg < 32; ++wg)
        for (int j = 0; j < 4; ++j)
            for (int kb = 0; kb < 8; ++kb)
                for (int L = 0; L < 32; ++L)
                    for (int r = 0; r < 4; ++r) {
                        const int rho = (L >> 2) + ((r & 1) ? 8 : 0);
                        const int y = (2 * j + rho / 8) * 256 + 8 * wg + (rho % 8);
                        const int kbase = ((r >> 1) ? 128 : 0) + 32 * (L & 3);
                        uint32_t word = 0;
                        for (int i = 0; i < 32; ++i) {
                            const int kap = kbase + i, W = 8 * kb + kap / 32, beta = kap % 32;
                            const int wgx = W / 2, hp = W % 2, lam = beta / 8, rk = beta % 8;
                            const int x = 256 * rk + 8 * wgx + 4 * hp + lam;   // byte = lane 8wg+4h'+lam, bit = round
                            uint32_t bit = 0;
                            if (x < sca::N && y < sca::N && x != y) bit = uint32_t((g.J[x][y >> 6] >> (y & 63)) & 1u);   // J[x][y], as run_trial
                            word |= bit << i;
                        }
                        A[((size_t(wg * 4 + j) * 8 + kb) * 32 + L) * 4 + r] = word;
                    }
    return A;
}

// Integer-range guard: the kernel evaluates z = s*h*65536 + q -/+ corr in int32. Prove no overflow for every step.
void check_ranges(const sca::Tables& tb) {
    const int64_t lim = (int64_t(1) << 31) - 1;
    for (size_t t = 0; t < tb.fourT.size(); ++t) {
        if (tb.fourT[t] <= 0) throw std::runtime_error("fourT must be positive");
        if (tb.kcorr[t] < 0) throw std::runtime_error("negative kcorr not supported");
        const int64_t corr_max = ((int64_t(sca::N) * tb.kcorr[t]) >> 8) + std::llabs(int64_t(tb.kconst[t])) + 1;
        const int64_t zmax = int64_t(sca::N - 1) * 65536 + std::llabs(int64_t(tb.q[t])) + corr_max + int64_t(tb.fourT[t]);
        if (zmax >= lim) throw std::runtime_error("table range exceeds int32 kernel arithmetic at step " + std::to_string(t));
    }
}

int64_t scalar_cut(const sca::Graph& g, const uint64_t* bits) {  // independent scorer (as the V80 host)
    int64_t c = 0;
    for (int x = 0; x < sca::N; ++x) {
        const bool sx = (bits[x >> 6] >> (x & 63)) & 1;
        for (int y = x + 1; y < sca::N; ++y)
            if (sx != (((bits[y >> 6] >> (y & 63)) & 1) != 0)) c += -sca::jval(g, x, y);
    }
    return c;
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
    int cs = 4, groups = 1;
    bool ons = false, lin = true;
    void* fn = nullptr;
    int threads = 0;
    size_t smem = 0;
    int max_clusters = 0;
    template <int CS, int G> void pick() {
        if (ons) fn = reinterpret_cast<void*>(kern::sca_kernel<CS, G, true, true>);
        else if (lin) fn = reinterpret_cast<void*>(kern::sca_kernel<CS, G, false, true>);
        else fn = reinterpret_cast<void*>(kern::sca_kernel<CS, G, false, false>);
        smem = kern::Cfg<CS, G>::TAB_OFF;
    }
    void setup(int S) {
        if (cs == 4 && groups == 1) pick<4, 1>();
        else if (cs == 4 && groups == 2) pick<4, 2>();
        else if (cs == 8 && groups == 1) pick<8, 1>();
        else if (cs == 8 && groups == 2) pick<8, 2>();
        else throw std::runtime_error("unsupported --cs/--groups (4|8, 1|2)");
        threads = 32 * (32 / cs);
        smem += size_t((12 * S + 15) / 16) * 16 + size_t(8) * S;
        CK(cudaFuncSetAttribute(fn, cudaFuncAttributeMaxDynamicSharedMemorySize, int(smem)));
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
        std::string graph = "data/K2000.bin", output, mode = "run";
        int chains = 1, batches = 1, cs = 4, groups = 1, trace = 0, clusters_arg = 0, warmup = 1;
        uint64_t seed = 20261004; uint32_t offset = 0;
        for (int a = 1; a < argc; ++a) {
            std::string k = argv[a];
            if (v80host::parse_schedule_flag(sc, k, a, argc, argv)) continue;
            if (k == "--info") { mode = "info"; continue; }
            auto val = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++a]); };
            if (k == "--graph") graph = val();
            else if (k == "--output") output = val();
            else if (k == "--chains") chains = std::stoi(val());
            else if (k == "--clusters") clusters_arg = std::stoi(val());
            else if (k == "--batches") batches = std::stoi(val());
            else if (k == "--cs") cs = std::stoi(val());
            else if (k == "--groups") groups = std::stoi(val());
            else if (k == "--trace") trace = std::stoi(val());
            else if (k == "--warmup") warmup = std::stoi(val());
            else if (k == "--seed") seed = std::stoull(val());
            else if (k == "--offset") offset = uint32_t(std::stoul(val()));
            else throw std::runtime_error("unknown argument " + k);
        }
        host::Device dev = host::open_device();
        const int S = sc.S;
        if (S < 2 || S > 4096) throw std::runtime_error("steps out of range");
        sca::Tables tb = v80host::build_tables(sc);
        host::check_ranges(tb);
        bool ons = false; for (int t = 0; t < S; ++t) ons |= tb.kcorr[t] != 0;
        host::Launcher ln; ln.cs = cs; ln.groups = groups; ln.ons = ons; ln.lin = ons || trace; ln.setup(S);
        const int per_cluster = 4 * groups;
        if (mode == "info") {
            std::cout << "{\"device\":\"" << dev.prop.name << "\",\"uuid\":\"" << dev.uuid << "\",\"sms\":" << dev.prop.multiProcessorCount
                      << ",\"onsager_kernel\":" << (ons ? "true" : "false") << ",\"nlin_kernel\":" << (ln.lin ? "true" : "false") << ",\"cs\":" << cs << ",\"groups\":" << groups << ",\"threads\":" << ln.threads << ",\"smem\":" << ln.smem
                      << ",\"max_active_clusters\":" << ln.max_clusters << ",\"chains_per_wave\":" << ln.max_clusters * per_cluster << "}\n";
            return 0;
        }
        if (output.empty()) throw std::runtime_error("--output PREFIX required");
        if (std::ifstream(output + ".trials.jsonl").good()) throw std::runtime_error("refusing to overwrite " + output);
        const int clusters = clusters_arg > 0 ? clusters_arg : (chains + per_cluster - 1) / per_cluster;
        if (chains < 1 || chains > clusters * per_cluster) throw std::runtime_error("chains must be in 1..clusters*4*groups");

        sca::Graph g = sca::load_graph(graph);
        if (g.sumw != -1040) throw std::runtime_error("unexpected graph (sumw != -1040)");
        std::vector<uint32_t> A = host::build_A(g);

        uint4* dA; int32_t *dFourT, *dQ, *dKconst, *dSh, *dTrace = nullptr; int64_t *dKcorr, *dFlips; uint8_t* dSpins;
        CK(cudaMalloc(&dA, A.size() * 4)); CK(cudaMemcpy(dA, A.data(), A.size() * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dFourT, S * 4)); CK(cudaMemcpy(dFourT, tb.fourT.data(), S * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dQ, S * 4)); CK(cudaMemcpy(dQ, tb.q.data(), S * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dKconst, S * 4)); CK(cudaMemcpy(dKconst, tb.kconst.data(), S * 4, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dKcorr, S * 8)); CK(cudaMemcpy(dKcorr, tb.kcorr.data(), S * 8, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&dSpins, size_t(chains) * 256)); CK(cudaMalloc(&dSh, size_t(chains) * 256 * 4)); CK(cudaMalloc(&dFlips, size_t(chains) * 8));
        if (trace) CK(cudaMalloc(&dTrace, size_t(chains) * S * 4));
        kern::Params p{};
        p.A = dA; p.fourT = dFourT; p.q = dQ; p.kcorr = dKcorr; p.kconst = dKconst; p.S = S; p.seed = seed; p.nchains = chains;
        p.do_trace = trace; p.out_spins = dSpins; p.out_sh = dSh; p.out_flips = dFlips; p.out_trace = dTrace;
#ifdef SCA_PROFILE
        long long* dProf; CK(cudaMalloc(&dProf, size_t(clusters) * cs * 8 * 8)); CK(cudaMemset(dProf, 0, size_t(clusters) * cs * 8 * 8)); p.prof = dProf;
#endif
        cudaStream_t st; CK(cudaStreamCreateWithFlags(&st, cudaStreamNonBlocking));
        cudaEvent_t e0, e1; CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));
        // warm-up launches on trial ids far outside every cohort (results discarded)
        for (int i = 0; i < warmup; ++i) { p.trial_offset = 0xF0000000u + uint32_t(i) * 65536u; ln.launch(clusters, p, st); }
        CK(cudaStreamSynchronize(st));

        std::ofstream trials(output + ".trials.jsonl"), bat(output + ".batches.jsonl"), spins(output + ".spins.bin", std::ios::binary);
        std::ofstream traces;
        if (trace) traces.open(output + ".trace.bin", std::ios::binary);
        if (!trials || !bat || !spins) throw std::runtime_error("cannot create output files");
        const std::string status_before = host::gpu_status();
        std::vector<uint8_t> hs(size_t(chains) * 256); std::vector<int32_t> hsh(size_t(chains) * 256), htr(trace ? size_t(chains) * S : 0);
        std::vector<int64_t> hfl(chains);
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
            CK(cudaMemcpy(hsh.data(), dSh, hsh.size() * 4, cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(hfl.data(), dFlips, hfl.size() * 8, cudaMemcpyDeviceToHost));
            if (trace) CK(cudaMemcpy(htr.data(), dTrace, htr.size() * 4, cudaMemcpyDeviceToHost));
            int bs = 0;
            for (int c = 0; c < chains; ++c) {
                uint64_t bits[32] = {0}; int64_t sh = 0;
                for (int l = 0; l < 256; ++l) {
                    const uint8_t v = hs[size_t(c) * 256 + l];
                    for (int k = 0; k < 8; ++k) {
                        const int i = k * 256 + l;
                        if (i < sca::N && ((v >> k) & 1)) bits[i >> 6] |= uint64_t(1) << (i & 63);
                    }
                    sh += hsh[size_t(c) * 256 + l];
                }
                const int64_t cut_dev = (g.sumw + sh / 2) / 2;
                const int64_t cut_host = host::scalar_cut(g, bits);
                if (cut_dev != cut_host) ++mism;
                const bool ok = cut_host >= 33000;
                succ += ok; bs += ok;
                const uint32_t tid = p.trial_offset + uint32_t(c);
                trials << "{\"batch\":" << b << ",\"chain\":" << c << ",\"trial_id\":" << tid << ",\"cut\":" << cut_host << ",\"device_cut\":" << cut_dev
                       << ",\"flips\":" << hfl[c] << ",\"success\":" << (ok ? "true" : "false") << "}\n";
                spins.write(reinterpret_cast<const char*>(bits), sizeof bits);
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
        const double Psec = 1.0 - std::pow(1.0 - pt, chains);
        std::ofstream sum(output + ".summary.json");
        sum << std::setprecision(10) << "{\"device\":\"" << dev.prop.name << "\",\"uuid\":\"" << dev.uuid << "\",\"cs\":" << cs << ",\"groups\":" << groups
            << ",\"clusters\":" << clusters << ",\"max_active_clusters\":" << ln.max_clusters << ",\"threads\":" << ln.threads << ",\"smem\":" << ln.smem
            << ",\"chains\":" << chains << ",\"batches\":" << batches << ",\"seed\":" << seed << ",\"trial_offset\":" << offset
            << ",\"onsager_kernel\":" << (ons ? "true" : "false") << ",\"nlin_kernel\":" << (ln.lin ? "true" : "false")
            << ",\"steps\":" << S << ",\"t0\":" << sc.t0 << ",\"t1\":" << sc.t1 << ",\"q\":" << sc.q << ",\"lambda\":" << sc.lam
            << ",\"ramp\":" << (sc.ramp ? "true" : "false") << ",\"tec_jv\":" << sc.jv << ",\"tecT_kappa\":" << sc.kappa << ",\"trace\":" << trace
            << ",\"successes\":" << succ << ",\"trials\":" << int64_t(chains) * batches << ",\"p\":" << pt << ",\"batches_with_success\":" << batch_succ
            << ",\"P_batch\":" << Pb << ",\"device_ms_mean\":" << tb_mean << ",\"device_ms_min\":" << dev_ms_min << ",\"device_ms_max\":" << dev_ms_max
            << ",\"host_ms_mean\":" << host_ms_sum / batches << ",\"tts99_primary_ms\":" << tts(Pb, tb_mean)
            << ",\"tts99_secondary_ms\":" << tts(Psec, tb_mean) << ",\"cut_mismatches\":" << mism
            << ",\"gpu_before\":\"" << status_before << "\",\"gpu_after\":\"" << status_after << "\"}\n";
        std::cout << "chains=" << chains << " batches=" << batches << " p=" << pt << " P_batch=" << Pb << " t_batch=" << tb_mean
                  << " ms  TTS_primary=" << tts(Pb, tb_mean) << " ms  TTS_secondary=" << tts(Psec, tb_mean) << " ms  cut_mismatch=" << mism << "\n";
        return mism ? 2 : 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
