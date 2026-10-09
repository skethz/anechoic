## ReAIM benchmark suite at the software level: all ten methods on Max-Cut, graph partitioning (GPP) and TSP

Folder: `research/reaim_benchmarks_20261008/`. Heavy runs were done on gpu-host under `/scratch/USER/anechoic_cpu_20261008` and copied back.

- The protocol was frozen at **01:48:15Z**, before any held-out run, and hashed in `PROTOCOL.sha256`. Every frozen hash still verifies.
- The per-problem results for all ten methods finished at **04:18 CEST**. The precision study, the export and the C++ bit-check finished at **04:39 CEST**, inside both of your deadlines.
- Everything below is software-model output, not hardware.

### Bottom line

- **The Onsager corrections help on Max-Cut but not on penalty-encoded problems.** On G1–G20, Onsager-κT needs **2.11× fewer steps than plain SCA**, and fewer on 20 of 20 instances. On GPP and TSP neither correction lifts the synchronous engine rules above plain SCA or TEC, and all four engine rules fail.
- **GPP: the four engine rules fail.**
  - Plain SCA, TEC, κT and online end in a balanced partition in only **3–6% of runs** at S = 4096. SA reaches 99.8% and bSB 99.3%.
  - The few balanced results are no better than a random bisection (quality 0.645 against a random baseline of 0.644).
  - **Why:** the dense balance penalty makes every spin flip together. Plain SCA flips 84% of spins per step on G1 and in the median run ends with all 800 spins on one side.
  - Retuning the penalty weight (P = 1, 2, 8; post hoc) does not fix this.
- **TSP: bSB is the only strong method.**
  - bSB reaches mean quality **0.976** (ARPD 2.5%) and reaches the 1% target on 3 of 6 instances. It also beats every value in ReAIM's Table VI.
  - The four engine rules always produce valid tours, but only at quality 0.475–0.482 (ARPD ≈ 113%). Random tours score 0.367.
  - SA, ReAIM, APC and dSB sit at 0.56–0.59.
- **Coupling bits needed (pre-registered study, bit-exact hardware arithmetic):**
  - Max-Cut needs K = 2 (it is exactly ternary).
  - TSP needs **K = 6** for all four rules (K ≤ 3 gives 0% valid tours; K = 4 gives 15–97%).
  - For GPP the study is uninformative, because the rules already fail at full precision.
- **Hardware export.** It is written: 45 instance/K pairs and 180 schedules. The C++ golden reference reproduces my model **bit for bit on 356 of 356 trials**.

### What ReAIM evaluated (checked in its text) and what I chose

- **ReAIM's setup.**
  - Instances: MCP G1–G20, GPP G1–G5 and G14–G17 (N = 800), TSP gr17, gr21, gr24, fri26, bayg29, bays29.
  - Iterations: 4096, 4096 and 8192.
  - Metrics: best of 20 runs for MCP and GPP; ARPD (mean over 20 runs) for TSP. Fig. 3 plots a "Normalized Cut Value" without saying how it is normalized.
  - ReAIM does **not** give the GPP/TSP formulas, the penalty weights, or how infeasible solutions are counted. For TSP it cites Ayodele (EvoCOP 2022), whose best TSP penalty is "MQC" (the maximum distance).
- **My pre-declared choices.**
  - **Formulations:** Lucas 2014 Ising forms with integer couplings (`lucas2014ising`).
  - **GPP penalty:** A = B = 1, giving couplings {−4, −3}. ReAIM does not state it; I inferred it as the only value consistent with ReAIM's 13-bit QUBO / 4-bit Ising widths.
  - **TSP:** one-hot with n² spins; penalty A = maximum distance.
  - **Field (bias) handling:** identical for all ten methods. The bias goes into the initial fields, equivalent to a frozen extra spin fixed at +1. The dense GPP penalty is carried exactly as a uniform coupling.
  - **Onsager-online scaling:** the coefficient is scaled by the mean row sum of J² over N (the Σ J² f′ rule). **Its limit:** with dense penalties J² is dominated by the penalty terms, and one global coefficient cannot represent the collective feedback that breaks the engines.
- **Data.** TSPLIB comes from the Heidelberg TSPLIB95 site (`reinelt1991tsplib`); two downloads were byte-identical, and SHA-256 values are in `data/tsplib_manifest.json`. All six optima were confirmed three ways: TSPLIB's list, the shipped optimal tours, and my own exact integer program. G-set data and hashes are the G-set study's (`helmberg2000spectral`, `ye_gset`).
- **Verification.**
  - Formulations checked by brute force over all states (GPP up to 14 nodes; TSP n = 3 and 4).
  - With zero bias, all kernels are bit-identical to the G-set study's `engine.py` / `others.py`. With bias, they are bit-identical to dense references.
  - Laptop and gpu-host are bit-identical on a 33-case fingerprint.
- **Procedure.**
  - Every method got 32 grid points, placed by an exploratory calibration on synthetic instances only.
  - Budgets: GPP S ∈ {256 … 4096}; TSP S ∈ {512 … 8192}.
  - 64-run pilots; selection by pilot mean quality, with infeasible runs scoring 0; 256-run held-out finals.
  - For Max-Cut, the four engine rules are reused unchanged from the G-set study. The six other methods were extended to the 15 missing instances with `run_gset.py`'s own grids, seeds and code.
- **GPP references.** These are my own SA-plus-swap values, computed before the held-out runs; they are not proven optima. Each is at or below ReAIM's Table V minimum: G1 7590, G2 7581, G3 7579, G4 7588, G5 7583, G14 1089, G15 1091, G16 1071, G17 1058. No method found a better cut.

### Compact paper table (8 rows)

Quality is ReAIM-style normalized quality (infeasible = 0) at ReAIM's budgets: MCP S = 4000, GPP 4096, TSP 8192. "Steps" is MCS99, the steps needed to reach the 1% target with 99% confidence: geometric mean of the per-instance minimum over S, over solved instances (solved/total).

| | SA | SCA | TEC | APC | ReAIM | aSB | bSB | dSB | Ons-κT | Ons-onl |
|---|---|---|---|---|---|---|---|---|---|---|
| Max-Cut G1–G20: quality | 0.997 | 0.991 | 0.993 | 0.994 | 0.976 | 0.984 | 0.989 | 0.995 | 0.994 | 0.991 |
| MCP: steps | 1,157 (20/20) | 8,581 (20/20) | 7,142 (20/20) | 4,040 (20/20) | 20,456 (17/20) | 3,276 (13/20) | 1,476 (17/20) | 2,309 (20/20) | 4,063 (20/20) | 5,140 (20/20) |
| GPP (9): quality | 0.985 | 0.039 | 0.034 | 0.616 | 0.856 | 0.108 | 0.980 | 0.023 | 0.031 | 0.017 |
| GPP: feasible | 1.00 | 0.06 | 0.06 | 0.91 | 0.88 | 0.11 | 0.99 | 0.03 | 0.06 | 0.03 |
| GPP: steps | 23,729 (8/9) | – | – | – | 182,720 (6/9) | 44,975 (5/9) | 12,833 (9/9) | – | – | – |
| TSP (6): quality | 0.567 | 0.476 | 0.482 | 0.560 | 0.586 | 0.000 | 0.976 | 0.563 | 0.477 | 0.475 |
| TSP: feasible | 1.00 | 1.00 | 1.00 | 1.00 | 0.92 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| TSP: steps | – | – | – | – | – | – | 21,621 (3/6) | – | – | – |

To read the quality scale, a random solution scores about 0.82 (MCP +1 graphs), 0.79 or 0.46 (GPP random or planar-like) and 0.34–0.45 (TSP); see `data/random_baselines.json`.

### Max-Cut G1–G20 (ten methods)

At S = 4000:

| Method | quality | best-of-20 | P(target) > 0 | solved | steps (geo-mean) | runs at best-known cut |
|---|---|---|---|---|---|---|
| SA | 0.9966 | 0.9996 | 20/20 | 20/20 | 1,157 | 763 |
| dSB | 0.9953 | 0.9991 | 20/20 | 20/20 | 2,309 | 320 |
| APC-SCA | 0.9942 | 0.9985 | 20/20 | 20/20 | 4,040 | 228 |
| Onsager-κT | 0.9941 | 0.9982 | 20/20 | 20/20 | 4,063 | 246 |
| Onsager-online | 0.9914 | 0.9975 | 20/20 | 20/20 | 5,140 | 205 |
| TEC | 0.9928 | 0.9979 | 20/20 | 20/20 | 7,142 | 152 |
| plain SCA | 0.9912 | 0.9975 | 20/20 | 20/20 | 8,581 | 142 |
| bSB | 0.9895 | 0.9918 | 15/20 | 17/20 | 1,476 | 47 |
| aSB | 0.9835 | 0.9872 | 11/20 | 13/20 | 3,276 | 2 |
| ReAIM ASA | 0.9760 | 0.9902 | 17/20 | 17/20 | 20,456 | 11 |

Paired step ratios over 20 instances:
- plain/κT 2.11× (κT fewer on 20/20);
- plain/online 1.67× (17/20);
- TEC/κT 1.76× (20/20);
- κT/SA 3.51× (SA fewer on 20/20);
- κT/dSB 1.76× (dSB fewer on 14/20).

The per-instance tables (20 × 10) are in `analysis.txt`. Steps are not time: an SA sweep is N sequential updates, an SCA step is one parallel update.

### GPP (G1–G5, G14–G17)

At S = 4096:

| Method | quality | feasible | quality when feasible | best-of-20 | solved | steps |
|---|---|---|---|---|---|---|
| SA | 0.985 | 0.998 | 0.987 | 0.994 | 8/9 | 23,729 |
| bSB | 0.980 | 0.993 | 0.987 | 0.992 | 9/9 | 12,833 |
| ReAIM ASA | 0.856 | 0.885 | 0.967 | 0.985 | 6/9 | 182,720 |
| APC-SCA | 0.616 | 0.915 | 0.662 | 0.686 | 0/9 | ∞ |
| aSB | 0.108 | 0.111 | 0.969 | 0.971 | 5/9 | 44,975 |
| plain SCA | 0.039 | 0.059 | 0.646 | 0.436 | 0/9 | ∞ |
| TEC | 0.034 | 0.059 | 0.645 | 0.425 | 0/9 | ∞ |
| Onsager-κT | 0.031 | 0.057 | 0.645 | 0.294 | 0/9 | ∞ |
| dSB | 0.023 | 0.027 | 0.840 | 0.337 | 0/9 | ∞ |
| Onsager-online | 0.017 | 0.030 | 0.644 | 0.152 | 0/9 | ∞ |

- **Per-instance quality** at S = 4096 (SA / bSB / ReAIM / APC / best engine rule):
  - G1–G5: 0.990–0.997 / 0.982–0.995 / 0.81–0.90 / 0.75–0.80 / ≤ 0.09.
  - G14–G17: 0.973–0.978 / 0.950–0.983 / 0.82–0.87 / 0.40–0.42 / ≤ 0.06.
- **Best balanced cut over all finals.**
  - SA and bSB: equal to R or within 1–12 edges of it.
  - Noise-free ReAIM: G1 7616, G2 7606, G3 7593, G4 7618, G5 7621, G14 1099, G15 1104, G16 1096, G17 1070. These are close to ReAIM's own Table V ASA values (7645 … 1075).
  - Engine rules: about 9,300–9,500 on G1–G5 and about 2,200 on G14–G17, roughly a random bisection.
- **GPP steps caveat.** G14–G17 are barely reached by SA (1 hit in 256 runs at S = 4096, so about 4.8 million steps). bSB reaches all nine.

### TSP (TSPLIB)

At S = 8192:

| Method | quality | feasible | best-of-20 | ARPD (%) over valid tours | at optimum |
|---|---|---|---|---|---|
| bSB | 0.976 | 1.00 | 0.986 | 2.5 | 15 runs |
| ReAIM ASA | 0.586 | 0.92 | 0.747 | 59.3 | 0 |
| SA | 0.567 | 1.00 | 0.667 | 79.1 | 0 |
| dSB | 0.563 | 1.00 | 0.672 | 80.0 | 0 |
| APC-SCA | 0.560 | 1.00 | 0.682 | 80.9 | 0 |
| TEC | 0.482 | 1.00 | 0.574 | ≈ 113 | 0 |
| Onsager-κT | 0.477 | 1.00 | 0.579 | 113.0 | 0 |
| plain SCA | 0.476 | 0.998 | 0.571 | 113.2 | 0 |
| Onsager-online | 0.475 | 0.999 | 0.565 | 113.7 | 0 |
| aSB | 0.000 | 0.00 | 0 | – | 0 |

Per-instance ARPD (%) at S = 8192, against ReAIM's Table VI:

| Instance | bSB | ReAIM (ours) | SA (ours) | SCA | κT | online | ReAIM paper ASA | ReAIM paper Neal | ReAIM paper best column |
|---|---|---|---|---|---|---|---|---|---|
| gr17 | 0.2 | 40.4 | 49.6 | 77.1 | 78.2 | 77.9 | 16.4 | 14.9 | 10.3 |
| gr21 | 6.2 | 58.9 | 78.9 | 113.8 | 111.3 | 111.4 | 32.1 | 33.0 | 10.7 |
| gr24 | 3.8 | 60.4 | 81.1 | 114.2 | 116.8 | 119.0 | 21.8 | 30.8 | 17.2 |
| fri26 | 2.0 | 67.4 | 90.0 | 124.4 | 125.5 | 124.0 | 44.0 | 29.4 | 29.4 |
| bayg29 | 1.4 | 63.2 | 90.4 | 121.4 | 119.8 | 122.5 | 31.0 | 26.3 | 26.3 |
| bays29 | 1.2 | 65.3 | 84.9 | 128.2 | 126.3 | 127.3 | 28.5 | 31.3 | 21.3 |

bSB reaches the 1% target on gr17 (P = 1.0 at S = 8192), gr21 and, barely, bays29. No other method reaches it on any instance.

### Hypotheses (pre-registered)

| ID | Result |
|---|---|
| R1 | **Pass.** Each engine rule is below SA on GPP on 9/9 instances |
| R2 | **Pass.** Neither Onsager form exceeds 50% feasibility on any GPP instance |
| R3 | **Pass.** bSB is best on TSP on 6/6 instances |
| R4 | **Pass.** Both Onsager forms are within 0.05 of plain SCA and TEC on 6/6 TSP instances |
| R5 | **Fail.** SA and dSB need fewer steps than every engine rule on 13 of 20 MCP instances (15 needed) |

### Coupling-precision study (pre-registered)

Setup:
- the four engine rules, run in a port of `sca_ref_bias.hpp` that is bit-exact to it;
- 1,024 paired trials per K, at each rule's selected schedule;
- tolerance: mean quality loss of at most 0.005 against full precision.

| Problem | Full precision (all four rules) | K = 2 | K = 3 | K = 4 | K = 6 / 8 | Smallest K within tolerance |
|---|---|---|---|---|---|---|
| MCP (20) | 0.987–0.990 | identical | ≈ full | ≈ full | ≈ full | **2** (all rules, all instances) |
| TSP (6) | 0.475–0.483, all valid | 0 valid | 0 valid | 0.33–0.36 (67–72% valid) | ≈ full | **6** (all four rules; 5/6 instances need 6, gr17 needs 4) |
| GPP (9) | 0.02–0.05 | graph information lost | noise | noise | noise | not informative (rules fail at full precision) |

- **TSP.** Rounding the distance couplings breaks the exact cancellation between the large bias and the coupling sums, which is why low K produces no valid tours.
- **GPP at K = 2.** The couplings {−4, −3} both round to −1, so the graph is lost.
- **Bias widths.**
  - Full-precision TSP biases need 17–18 signed bits, more than the spec's B_BITS = 16. At K = 8 they need 14–15 bits.
  - Full-precision TSP couplings span 119–257 distinct levels with max |J| = 560–1,730. At K = 8, 50–58 levels survive.
- **Model cross-check.** Float model against hardware arithmetic at full precision, mean difference: MCP −0.0001, TSP +0.0004, GPP −0.011 (GPP qualities are near 0 and the selection was noisy).

### Hardware export (`hw_export/`)

- **Contents.** For each GPP and TSP instance at K = 2, 4, 8:
  - `.jint8` couplings (SCAJINT8 format) and `.bias.bin` (int32);
  - `selected_configs.json` with, per rule: q, T0, T1 (final T), S, λ_eff or κ_eff, ramp; the range checks; the objective mapping and decoding; the measured precision-study result; SHA-256 of every file.
- **C++ check.** `sca_ref_bias.hpp` built with c++ −O2 −std=c++17 −DSCA_LANES=256 on gpu-host: **356/356 trials bit-identical** (spins, flips, the per-step count of spins in the linear band, Σ s·h, Σ b·s). Two schedules were skipped because their tables overflow int32: Onsager-online at K = 8 on G3 and G17.
- **Range limits the hardware teams should know about.**
  - 20 of 24 TSP K = 8 schedules exceed the FPGA's 27-bit 4T lane constant (T up to about 1,016, limit 512).
  - So do 7 of 36 GPP K = 8 schedules.
  - The |q| + |corr| < 2^29 bound fails for 1 of 36 GPP schedules at K = 2 and 3 of 36 at K = 4.

### Exploratory, post hoc (labelled; hashed before running)

- **ReAIM with F = min** (ReAIM's own Step-4 option; my implementation, inherited from the G-set study, fixes F = max).
  - GPP improves to 0.967–0.985 on G1–G5 and 0.905–0.922 on G14–G17 (from 0.81–0.90 and 0.82–0.87).
  - TSP gets worse (0.35–0.51 against 0.55–0.69), so F = max stays the right choice there, despite ReAIM's remark that TSP leans towards min.
- **GPP penalty P ∈ {1, 2, 8}** (on G1 and G14).
  - The engine rules stay at quality ≤ 0.21 and feasibility ≤ 0.25 at P = 1 and 2.
  - SA is best at P = 2 (G1 0.999); bSB is best at P = 8 (G1 0.998).
- **Bits needed by the methods that do solve GPP/TSP** (SA and bSB, float model).
  - GPP: K ≥ 3 is enough for SA.
  - TSP: SA needs K ≥ 6, and bSB still loses 1–7% at K = 8 on 4 of 6 instances (gr17 0.969 against 0.998 at full precision). Near-optimal tours need more than 8 bits under max-scaling.

### Caveats

1. **Software model only, no timing.** Today's V80 stores ±1 couplings. GPP and TSP need the multi-bit engine, and TSP also needs the bias path.
2. **Penalties are fixed by rule** (GPP A = B = 1, inferred; TSP A = maximum distance). Other penalties change the results; the engine-rule failure on GPP held at P = 1, 2 and 8.
3. **ReAIM is noise-free with F = max**, with T starting at 1 rather than ReAIM's TSP 0.5 → 0.1. With F = min (post hoc), its GPP results are much better.
4. **TSP gap with ReAIM's paper.** My SA (ARPD about 79%) and ReAIM (about 59%) are much worse than ReAIM's reported Neal (15–33%) and ASA (16–44%). I could not resolve this. Possible reasons: Ayodele's encoding used (n−1)² variables with city 1 fixed (I used n²), and the run and reads conventions are not stated.
5. **Grid edges** (`edge_report.txt`).
   - TSP engine rules select the largest T_fin in 28–30 of 30 selections, and the smallest T0 in 17–25.
   - APC on GPP selects the largest q_reset in 30/45 and the largest T0 in 38/45.
   - ReAIM selects the smallest k set on TSP in 29/30.
   - Calibration suggested no large hidden gains, but this biases against those methods.
6. **GPP references and targets** come from my own heuristic, not proven optima.
7. **One log file is messy.** `logs/precision_part3.log` has null padding, because a duplicate runner was briefly started and killed. Results are unaffected: the seeds are deterministic, writes are atomic, and all 840 files hold 1,024 trials.
8. **No writes outside scratch or this folder.** Nothing went under `/home/USER` (the recent `~/.ccache` changes are not from my plain g++ build). The old Snowball path was not created. No hashed or other research folder was edited.

### Recommended paragraph

> **Beyond Max-Cut (software model).** Following ReAIM's benchmark suite [10609617], we ran all ten methods on balanced graph partitioning (G-set G1–G5, G14–G17 [helmberg2000spectral; ye_gset]) and TSP (six TSPLIB instances [reinelt1991tsplib]) in Lucas's Ising forms [lucas2014ising] with fixed penalties, using the same pre-registered tuning (32-point grids, 64-run pilots, 256-run held-out finals). On Max-Cut G1–G20, the κ·T correction reached 99% of the best-known cut in 2.1× fewer steps than plain SCA, fewer on all 20 instances. On the penalty-encoded problems, all four synchronous SCA rules, with or without the correction, failed. On partitioning, the dense balance penalty makes every spin flip together, so only 3–6% of runs ended balanced, against ≥ 99% for SA and simulated bifurcation (bSB). On TSP they returned valid tours of mean normalized quality 0.48 (random tours: 0.37), against 0.98 for bSB. Such problems therefore need constraint-aware or sequential updates, not only more coupling precision; TSP also needs at least 6-bit couplings in the engine's arithmetic.

### Files

All in `research/reaim_benchmarks_20261008/`:

- **Protocol:** `PROTOCOL.md`, `PROTOCOL_ADDENDUM1.md`, `PROTOCOL.sha256`
- **Data:** `data/` (`tsplib_manifest.json`, `gpp_reference.json`, `instance_stats.json`, `random_baselines.json`, `tsplib/`), `sources/` (TSPLIB pages, Ayodele 2022)
- **Code:** `problems.py`, `solvers.py`, `grids.py`, `run_bench.py`, `hwmodel.py`, `hw_check.cpp`, `precision.py`, `export_hw.py`, `analyze.py`, `precision_analyze.py`, `explore_*.py`, `calib*.py`
- **Verification:** `verify_hope.json`, `verify_hw.json`, `verify_hw_frozen.json`, `repro_{laptop,hope}.json`, `hw_export/export_check.json`
- **Results:** `results/{gpp,tsp,mcp_ext}/`, `results_summary.json`, `analysis.txt`, `compact_table.md`, `edge_report.txt`, `precision/` (`schedules.json`, `runs/`, `precision_summary.json`, `precision_analysis.txt`), `explore/` (`explore_summary.txt`), `hw_export/`, `logs/`, `calib/`

Sources: [Ayodele 2022, arXiv:2206.11040](https://arxiv.org/abs/2206.11040); [TSPLIB95 (Heidelberg)](https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/).
