# Addendum 1: implementation files written after the freeze (no change to the frozen design)

Written on 8 October 2026 after `PROTOCOL.md` was frozen (01:48:15Z) and while the held-out runs (started 01:49:26Z) were running, **before any precision-study run and before any analysis output was seen**. The design is unchanged.

The following files implement the measures, the precision study and the export exactly as specified in `PROTOCOL.md`. Their SHA-256 values are appended to `PROTOCOL.sha256` before they are run.

- **Analysis of the held-out finals.**
  - `analyze.py` (all ten methods, MCP/GPP/TSP).
  - `precision_analyze.py` (the precision study; tolerance −0.005, K_min as pre-declared).
- **Precision study.** `precision.py`:
  - schedule selection (S_exp rule);
  - the quantized-schedule mapping;
  - 1,024 paired trials per K;
  - the seeds of `PROTOCOL.md`.
- **Export and checks.**
  - `export_hw.py`: the SCAJINT8 and int32 bias files, `selected_configs.json`, and the bit-for-bit check of the exported files against the C++ reference.
  - `fetch_v2.sh`: copies results back without overwriting frozen local files.
- **Re-check.** `verify_hw.py` was re-run on gpu-host with the frozen `hwmodel.py` (sha `bfc8ef66…`). Its output `verify_hw_frozen.json` is byte-identical to the frozen `verify_hw.json` (60/60 identical).

## Exploratory, post hoc (labelled as such in the report; not part of the pre-registered analysis)

- **`baselines.py`.** The expected normalized quality of a uniformly random solution (random spins, random exact bisection, random tour), used only to read the normalized scales (`data/random_baselines.json`).
- **`explore_penalty.py`.** GPP penalty sensitivity, P ∈ {1, 2, 8} against the protocol's P = 4, on G1 and G14 at S = 4096.
  - Same grids, pilot, selection and final rules as the protocol, with fresh seeds.
  - This was prompted by the calibration finding that synchronous updates fail on the dense balance penalty.
