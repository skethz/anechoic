# V80 Onsager-SCA engine: design spec (draft for review)

3 October 2026. Target: AMD Alveo V80 (`xcv80-lsva4737-2MHP-e-S`) on fpga-host.
- Reuses the `fpga/v80_snowball` AVED shell, host access path, timer convention and build/programming scripts.
- Algorithm: SCA (STATICA Eq. 7, uniform thresholds) with the Onsager echo correction and optional end ramp ([REPORT_ABCD.md](../../research/theory_ideas_20261003/REPORT_ABCD.md)).
- Goal: measure real flips per cycle, clock frequency and TTS, then replicate engines.

## Algorithm, per trial

1. Seed per-lane RNGs from (seed, trial). Draw random initial spins s.
2. Initial fields `h = J s`, by streaming all N rows.
3. For each step t = 0..S−1, with host-supplied tables `fourT[t]` (Q8.16), `q[t]` (Q16.16) and `kcorr[t] = λ_t / (2 T_{t−1})` (Q8.24):
   - **Decision.** For each spin i:
     - `corr = (n_lin(t−1) · kcorr[t]) >> 8` (Q16.16).
     - `z = s_i·h_i·2^16 + q[t] − corr·(s_i == sprev_i ? 1 : −1)`.
     - `r = ((u16 − 2^15) · fourT[t]) >> 16`, giving a uniform value in [−2T, 2T) in Q16.16.
     - Flip iff `z < r`. This is exactly Eq. 7, `clamp(z/4T + ½) < u`, up to the 16-bit quantisation of u.
     - `n_lin(t) = #{|z| < 2T}`.
   - **Commit.** sprev ← s, then flip the selected spins.
   - **Field update.** For each flipped spin x, add `2·J_xy·s_x(new)` to all h_y.
4. Score: `cut = (Σw + ½ Σ_i s_i h_i)/2`. Write the final spins, cut, total flips and cycle count to HBM.

There are no floating-point operations in the engine. The schedule (including λ ramps and q schedules) lives entirely in the host tables, so one bitstream covers plain SCA (kcorr = 0), Onsager SCA and the TEC variant (a constant-correction mode bit).

## Per-engine datapath (STATICA-style, flip-driven)

| Block | Choice | Rationale |
|---|---|---|
| Coupling store | 2,048 rows × 2,048 bits in 32 dual-port BRAM banks of 64 bits | One full row per port per cycle, so **2 flipped spins per cycle**. About 128 BRAM36 per engine. |
| Fields | 2,048 × 13-bit registers, each with a 3-input adder (two rows) | All fields updated every cycle, as in STATICA's LAUs |
| Flip list (DDSS) | Flip bitmap plus a word-occupancy bitmap; two find-first-set extractions per cycle | Matches STATICA's priority-circuit DDSS |
| Decision lanes | L = 128 initially (2,048/L = 16 cycles per step); one 16×24 DSP multiply and an xoshiro128** RNG per lane | L is a parameter; 256 lanes means 8 cycles per step |
| Counters | 64-bit cycle counter; flips and n_lin per step | Timing and model validation |

**Expected cycles per step:** `2048/L + ⌈F/2⌉ + pipeline` (roughly 20 cycles of overhead at L=128). The STATICA-fitted model is about 0.43 cycles per flip with ~0 per-step overhead, so with L=128 we expect it to be **slower at low flip counts**, which shows up near the end of the anneal. Measuring this is one of the first goals.

## Verification

1. A **C++ bit-exact reference** (`sca_ref.hpp`) uses the same fixed-point formats, RNG, lane order and tables.
   - Its quality must match the floating-point Python model within confidence intervals at the selected configurations. This guards against the STATICA "folded RNG" failure mode.
2. **HLS C simulation** must match the reference bit for bit: final spins, cut, flips and the per-step n_lin trace, for short and full schedules on K2000 and on small random graphs.
3. **RTL co-simulation** on a short schedule.
4. **On board:** UUID and timer checks, an HBM readback, and exact replay of seeded trials against the reference before any timed cohort.

## Milestones

| M | Deliverable |
|---|---|
| M1 | Reference model and quality check (laptop) |
| M2 | HLS kernel; C simulation equals the reference; synthesis estimates. On fpga-host under `/scratch/USER/sca_v80_20261003`. |
| M3 | Board integration, routed timing at 300 MHz (or best), programmed image |
| M4 | Physical cohort with measured cycles per flip and TTS against the cycle model and GPU |
| M5 | Replicate to 4/8/16 engines as resources and timing allow |

## Open choices for review

- Lane count L (128 vs 256): decision latency against DSP cost when replicating.
- 16-bit uniforms (STATICA also used 16-bit) vs 32-bit.
- Whether to also include the logistic kernel needed for exact-weight population annealing (C). It would cost a LUT-based logistic function in each lane.

## Implementation log

| Version | Change | HLS result (xcv80, estimates use the versal-aicore characterization fallback) |
|---|---|---|
| M1 | Bit-exact C++ reference | Matches the float model: plain 0.812 vs 0.824, Onsager 0.865 vs 0.866, ramp 0.726 vs 0.699. The first attempt exposed a threshold-scaling bug (`>>15` gave a [−4T, 4T) band); fixed to `>>16`. |
| v1 | Direct implementation | All loops II=1, but the UPDATE recurrence (word select → ctz → clear → ctz) was 8.5 ns: **Fmax ≈ 117 MHz**, 548k LUT. |
| v2 | Two flip extractors (port A: spins 0–1023, port B: 1024–2047) with next-word prefetch; diagonal handled by the invariant `h_stored = h_true − s` (decide adds +1, score adds N); one shared LANE loop (init / decide / score) and one ACCUM loop (init fields / update) | **Fmax ≈ 231 MHz**, 381k LUT (14%). New critical paths: the cross-lane score adder chain in LANE (4.3 ns) and select-chain `ctz` loops in ACCUM (3.8 ns). |
| v3 | `ctz` via isolate-lowest-bit plus OR-tree encoders; spin bit via `(word & onehot) != 0`; per-lane accumulators for the score and n_lin, reduced outside the II=1 loop | LANE fixed (4.3 → 2.4 ns). ACCUM 4.5 ns: the next-word prefetch was still in the cur recurrence, and HLS built the encoder as a select chain. **Fmax ≈ 221 MHz**, 397k LUT. |
| v4 | Encoder bits computed as independent OR reductions; addresses composed as `(word << 6) \| bit` instead of added; prefetched next-word register, refilled from occ only when consumed | HLS estimate **276 MHz**. RTL co-simulation PASS. **Vivado routed at 300 MHz** (WNS 0.000 ns, 0 failing endpoints) and at 250 MHz (+0.030 ns). Whole image 200k LUT (7.8%). On the board: 24/24 exact; results in [RESULTS_HW.md](RESULTS_HW.md). |
| v5 | **256 lanes** (8 rounds per step), **4 interleaved flip extractors** with two coupling copies (4 read ports), burst-loaded tables. The lane count is a compile-time macro (`SCA_LANES`, default 128 = v4); v5 uses 256. Statistically neutral (O1 p 0.932 vs 0.933 in the reference). Built as `SCA_VARIANT=v5` with separate HLS project, board tag, host and manifest. | Native bit-exact at 256 and 128 lanes; C simulation passes. HLS: II=1 everywhere, **Fmax ≈ 206 MHz**; 587k LUT (22%) and 563 BRAM (HLS estimates). The critical path is the ACCUM refill: occ → lowest bit → encoder → 8:1 word mux → 4 select levels (4.86 ns). |
| v5.1 | v5 + **queue-based refill**: at step start each extractor compacts its non-empty words into an ordered queue, so a refill is one indexed read. x and spin are computed unconditionally (no early return), removing select chains on the coupling address. Variant `v51`, source `sca_kernel_v51.cpp`. | HLS **297 MHz**, co-simulation PASS. **Vivado: 300 MHz WNS −0.269 ns, 275 MHz −0.060 ns (failed).** The worst paths are about 80% routing: the per-flip ±d broadcast to the 2,048 field adders, and the w → word → spin → d path. v5 at 250 MHz closed (+0.001 ns) and was measured: 0.275 cycles/flip, 31.1 cycles/step. |
| v5.2 | v5.1 + **one packed table array** (4 × 64-bit words per step, single burst stream; v5 missed bursts on the three 32-bit tables). New interface (`tab_in`); host `host_sca_multi.cpp` (1–5 engines). | HLS 297 MHz. **4-engine build at 250 MHz did not close**: post-placement WNS −1.8 ns on engine-internal LANE/ACCUM paths (stopped; 225 MHz fallback not attempted). The first attempt also exposed a packaging bug: the timer registers no longer fit the 7-bit register map. They now use undecoded `*_CTRL` placeholder words. |
| v5.3 | v5.2 + **spin-word queue** (the spin sign comes from a per-word register, with no dynamic word mux in the loop) + **registered per-field partial sums** (each iteration's 4-term coupling sum is added to h in the next iteration, with a drain add after the loop). Reorders additions only. | Native bit-exact; HLS **331 MHz** (ACCUM 3.02 ns); co-simulation PASS. Single engine at 300 MHz: post-placement WNS −0.40 ns (275 MHz retry queued). **2 engines at 250 MHz: post-placement WNS −1.68 ns (stopped).** Final single-engine results: 300 MHz WNS −0.176 ns, 275 MHz −0.638 ns. Both failed; run-to-run placement variation is about 0.5 ns. **3 engines with one pblock per SLR at 250 MHz: final WNS −0.137 ns** (near miss; the floorplan removed most of the cross-die penalty). Remaining violators are engine ACCUM/LANE paths and the shared control/memory SmartConnect crossings. **225 MHz: routed −0.014 ns, then a post-route phys_opt (AggressiveExplore) closed it to 0.000 ns. Image built, programmed and measured: best TTS99 0.211 ms (RESULTS_HW.md).** The same post-route pass on the 250 MHz build only reached −0.124 ns (TNS −6,483 ns) and was refused by the timing gate. Going faster needs a structural change: a deeper, replicated broadcast pipeline for the per-flip address and ±d fan-out. |

**Replication finding.** One v5-class engine closes 250 MHz with no margin (+0.001 ns). Two or four copies fail by 1.7–1.8 ns on engine-internal ACCUM/LANE paths.
- Each engine broadcasts its 4 extracted addresses and ±d values to 2 × 32 coupling BRAM banks and 2,048 field adders. Those resources are spread over the die's BRAM columns, so a second engine forces both to spread further.
- Fixing this needs per-engine floorplanning (pblocks) and/or a deeper, replicated broadcast pipeline: registered address/d copies per bank group, adding about 2 cycles per step. More copies of the current engine alone will not close.

Every version matched the reference bit for bit on the native testbench: 10 trials covering plain, Onsager with ramp and TEC, including a full 560-step schedule. v1 and v2 also passed Vitis HLS C simulation.

## v6: hand-written RTL engine (6 October 2026)

**Why it was needed.**
- The HLS engines are limited by HLS itself, not by the algorithm.
- Every pipelined loop (SEED, LANE, ACCUM) is outlined into its own module, and the register-resident state (fields, RNG states, spin words) is passed in and out by value. So the state exists 2–3 times: v5.3 has 356k FF per engine against about 100k bits of real state.
- The loop-exit handshakes drive 32k–43k-fanout clock enables through BUFG_FABRIC. These are the worst paths at 250–300 MHz: AXI write response → RNG CE; loop-exit compare → field CE.
- The 3-engine 250 MHz build had 137,640 failing endpoints, almost all CE pins.

**Structure (`src/v6/sca_core.v`).**
- **Groups.** There are 32 groups. Group b owns:
  - spins and fields y = b + 32j;
  - lanes l = b + 32m (lane l takes spins l + 256k in round k, i.e. field j = 8k + m);
  - coupling bank b (bit j of address x = J[x][b + 32j]) in two copies.

  The field update, the lane field multiplexers and the coupling reads are therefore all local to a group.
- **Extractors.** Four extractors drain one 64-bit word per round: bit 8m + q of extractor e comes from group e + 4q, lane m. Extraction starts while the remaining decision rounds are still running; flipped spin x = {k, p, e}.
- **Field update.** Field j adds m·(2·popcount(J ? P : Q) − nvalid), with registered partial sums. All state exists once.
- **Seeding and I/O.**
  - Lanes are seeded through a 256-stage chain from one pipelined splitmix64 unit, so there is no broadcast.
  - The interface is a stream: a 128-bit input carries the header, couplings and tables; a 64-bit output carries results and traces.
  - The HLS front end `src/v6/sca_mover.cpp` keeps the v5.2 register map and timers. It uses 2 memory masters per engine instead of 3.

**Cycle accounting.**
- The per-step overhead is max(about 23, about 13 + max-extractor flips) cycles, against 31 + 0.2754·flips in v5.
- Fitted on RTL simulation of the six protocol trials: cycles/trial = 0.2714·flips + 16.91·S + 1148, which is 0.77–0.82× of the v5 model.
- DRAIN = 2 is the exact minimum: simulation passes at 2 and fails at 1.

**v6.0 → v6.1 (timing only; cycles unchanged).**
- **v6.0.** Standalone at 300 MHz it reached WNS −0.78 ns. Failing endpoints:
  - 3,403: RNG-enable fanout;
  - 6,939: per-step constants into 33-bit lane adders;
  - about 1,100: the BRAM address mux;
  - 426: the score adder chain.
- **v6.1 fixes.**
  - registered and replicated enables and constants;
  - an exact 16-bit high/low split of z (two's-complement lexicographic compare against thr, U = twoT − 1 and L = 1 − twoT, with precomputed low-half comparisons);
  - a registered address mux;
  - a pipelined score path;
  - an input skid buffer.
- **v6.1 result.** +0.323 ns at 300 MHz, 0 of 272k endpoints failing. 128,118 LUT, 126,373 FF, 258 BRAM tiles and 282 DSP per engine.
