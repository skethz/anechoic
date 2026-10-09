# K = 4 + bias image (Amendment 4): build and board status

This file is updated as the watch proceeds. All times are UTC; KST is UTC + 9 h. fpga-host identity check: `fpga-host USER`, 15:32 UTC on 8 Oct.

## 250 MHz build (`board_mb4b_e6_250mhz`): failed in routing

- **Started** 07:08 UTC; pipeline exit 1 at 12:06 UTC.
- **Post-placement** physopt WNS +0.001 ns.
- **Routing.** Estimated congestion level 5 (32×32), with congested clusters at level 4. `route_design` failed with:
  - 4,679 signals unroutable because of congestion;
  - 4,054 node overlaps;
  - "Design is not legally routed".

## 225 MHz rebuild (`board_mb4b_e6_225mhz`): pre-declared fallback, started 12:06 UTC

**State at 15:32 UTC: still running, in routing (Phase 5, rip-up and reroute).**

- **Post-placement utilisation:**

| Resource | Used | Device |
|---|---:|---:|
| LUT | 2,209,779 | 85.8% |
| FF | 1,370,759 | 26.6% |
| Block RAM tiles | 3,183 | 85.1% |
| URAM | 1,566 | 81.4% |

- **Timing so far.**
  - Post-placement physopt WNS +0.200 ns.
  - Initial routing estimates congestion level 5 (32×32) again.
  - Intermediate routed timing: WNS −0.029 ns (TNS −0.704 ns), WHS +0.010 ns.
- **Board.** v6.4 headline image programmed and idle. The Amendment-4 chain is waiting for this image; it will not program after 02:00 UTC on 9 Oct.

## Updates

(appended below by the watch)

### 17:04 UTC: the 225 MHz rebuild routed and passed timing

- **Pipeline** exit 0 at 17:01:58 UTC.
- **Post-route timing:** setup **+0.002 ns**, hold +0.010 ns, 0 failing endpoints; signed off as "PASS for constrained synchronous timing".
  - The router went through −0.049 ns and −0.029 ns before it converged.
- **Routed utilisation:**

| Resource | Used | Device |
|---|---:|---:|
| LUT | 2,213,218 | 86.0% |
| FF | 1,370,759 | 26.6% |
| RAMB36 | 3,180 (+6 RAMB18; 3,183 Block RAM tiles) | 85.1% |
| URAM | 1,566 | 81.4% |
| DSP58 | 1,692 | 15.6% |

- **Per-SLR Block RAM:** 1,061 tiles each, which is 79.1% / 88.4% / 88.4%.
- **Image:** `board_mb4b_e6_225mhz`, logic UUID 015c977e7098d29dfc7ea51a2f37bf6c, kernel clock 225 MHz (224.999999775 nominal), PDI SHA-256 2dd9b530fc904c5d….
- **Amendment-4 chain:** took this image at 17:02:03 UTC. It backed up the staged v6.4 PDI (SHA-256 67e591f3…, `results/staged_image_backup_20261008T170203Z`) and started programming at 17:02:04 UTC.

### 18:55 UTC: Amendment 4 completed, all verdicts PASS, v6.4 restored

**Chain timeline (UTC):** programming 17:02–17:52; amendment 17:52:25–17:52:51; power 17:52–18:04; v6.4 restore 18:04–18:49 (UUID 8669ad38… matched, BAR2 OK); analysis 18:49. Programming finished before the 02:00 UTC cutoff.

**Records:** `results/multibit_20261007/board/board_mb4b_e6_225mhz/`, including `hwmb_board_mb4b_e6_225mhz_a4/analysis.txt`, `a4_summary.json`, the power directory, `build/`, `build_250mhz_failed/` and the logs.

| ID | Result | Verdict |
|---|---|---|
| A4-E | 19 sets, 156/156 trials bit-exact against `run_trial_bias` (spins, score word, flips, n_lin trace), all 6 engines, padding bits 0 | **PASS** |
| A4-R | K2000 with bias mode off: X5 and O4 each 2052/2052 trial ids with the same cut and flips as the 12 × v6.4 run; W1–W4 identical | **PASS** |
| A4-C | Launch cycles against v6.4: X5 0.99947, O4 0.99953 | **PASS** |
| A4-F | TSP valid tours among trial ids 0–1023 equal the export exactly: gr17 993, gr21 595, gr24 974, fri26 844, bayg29 158, bays29 678 | **PASS** |

**A4-E breakdown:**
- **Bias mode off:** 42/42 (K2000 W1–W4).
- **Bias mode on:** 114/114, made up of:
  - K2000 with b = 0 (W1b0, W4b0): 24/24;
  - random K = 4 matrices with random biases at n = 777 and 2048: 36/36;
  - GPP exports G1, G14, G17 (dense K = 4, n = 800): 18/18;
  - TSP exports with bias (n = 289–841, S = 8192): 36/36.

**K2000 TTS99 on this image** (6 engines at 225 MHz, primary): X5 0.0981 ms (p 0.338, t_R 0.0557 ms) and O4 0.1171 ms. The 12 × v6.4 image at 250 MHz gives 0.0503 and 0.0528 ms.

**Board power** (hwmon4; 150 s loads, 512 trials per launch, devices 99.7–99.8% busy):

| Load | Board W | Dynamic W | Notes |
|---|---:|---:|---|
| Idle | 67.3 (first); 68.6–70.5 between loads | | |
| X5 | 93.5 | 24.7 | 9.17 mJ to solution; 854 µJ per trial |
| bays29 TSP (n = 841, S = 8192) | 78.2 | 8.4 | |
| G1 GPP (all spins flipping) | 77.0 | 8.1 | |

**Sentence for §5.5:** With K = 4 and per-spin bias, the URAM read ports allow six engines; they fail routing at 250 MHz but close timing at 225 MHz (86% LUT, 81% URAM). That image is bit-exact against the golden model on all 156 checked trials (114 with bias), so the device runs six K = 4 engines at 225 MHz, while K = 8 would fit only three.
