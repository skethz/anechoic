# Erratum: `hw_export/selected_configs.json`, field `precision_study` (GPP cells)

Found 2026-10-08 ~05:10 CEST following the GPU study's report (gpu/multibit_bias_20261008): its K = 4 GPP runs on the export's own seeds did not reproduce 16 exported `precision_study` entries.

**Cause.** The GPP instances (G1-G5, G14-G17) share names with the Max-Cut instances G1-G20. In 49 GPP cells, the export's `precision_study` field holds the values of the **Max-Cut** cell with the same instance name, rule and K (an instance-name collision in `export_hw.py`'s lookup). Every mismatching cell equals the corresponding `mcp/<inst>/<rule>` entry of `precision/precision_summary.json` exactly.

**Scope.**
- Only this metadata field is affected. The matrices, biases, schedules, seeds and range checks in the export are unaffected; the GPU study rebuilt the K = 4 matrices and biases independently and they match exactly.
- `precision/precision_summary.json`, `precision/precision_analysis.txt` and all benchmark results are correct; the GPP cells there show failure at every K, consistent with the C++ golden reference.
- The manuscript uses only `precision_summary.json` / `precision_analysis.txt` and the benchmark finals, so no manuscript number changes.

**The export file is not edited** (hashes preserved). Use the values below, taken from `precision/precision_summary.json`, as the expected GPP results.

| GPP instance | K | rule | exported (wrong): quality / p_feasible / p_target | correct (precision_summary): quality / p_feasible / p_target |
|---|---|---|---|---|
| G1 | 2 | SCA | 0.9920 / 1.000 / 0.996 | 0.0177 / 0.022 / 0.000 |
| G1 | 2 | Onsager-kT | 0.9931 / 1.000 / 0.946 | 0.0270 / 0.034 / 0.000 |
| G1 | 2 | Onsager-online | 0.9861 / 1.000 / 0.984 | 0.0271 / 0.034 / 0.000 |
| G1 | 4 | SCA | 0.9919 / 1.000 / 0.996 | 0.0613 / 0.077 / 0.000 |
| G1 | 4 | Onsager-kT | 0.9931 / 1.000 / 0.946 | 0.0301 / 0.038 / 0.000 |
| G1 | 8 | SCA | 0.9919 / 1.000 / 0.996 | 0.0627 / 0.079 / 0.000 |
| G1 | 8 | Onsager-kT | 0.9931 / 1.000 / 0.946 | 0.0286 / 0.036 / 0.000 |
| G2 | 2 | TEC | 0.9941 / 1.000 / 0.995 | 0.0255 / 0.032 / 0.000 |
| G2 | 2 | Onsager-online | 0.9839 / 1.000 / 0.801 | 0.0201 / 0.025 / 0.000 |
| G2 | 4 | TEC | 0.9941 / 1.000 / 0.995 | 0.0255 / 0.032 / 0.000 |
| G2 | 8 | TEC | 0.9941 / 1.000 / 0.995 | 0.0263 / 0.033 / 0.000 |
| G3 | 2 | SCA | 0.9924 / 1.000 / 0.988 | 0.0262 / 0.033 / 0.000 |
| G3 | 2 | Onsager-kT | 0.9935 / 1.000 / 0.966 | 0.0254 / 0.032 / 0.000 |
| G3 | 4 | SCA | 0.9924 / 1.000 / 0.987 | 0.0698 / 0.088 / 0.000 |
| G3 | 4 | Onsager-kT | 0.9935 / 1.000 / 0.966 | 0.0285 / 0.036 / 0.000 |
| G3 | 4 | Onsager-online | 0.9847 / 1.000 / 0.979 | 0.0262 / 0.033 / 0.000 |
| G3 | 8 | SCA | 0.9924 / 1.000 / 0.987 | 0.0796 / 0.101 / 0.000 |
| G3 | 8 | Onsager-kT | 0.9935 / 1.000 / 0.966 | 0.0285 / 0.036 / 0.000 |
| G3 | 8 | Onsager-online | 0.9847 / 1.000 / 0.979 | 0.0309 / 0.039 / 0.000 |
| G4 | 2 | SCA | 0.9861 / 1.000 / 0.981 | 0.0208 / 0.026 / 0.000 |
| G4 | 2 | Onsager-kT | 0.9936 / 1.000 / 0.948 | 0.0232 / 0.029 / 0.000 |
| G4 | 4 | SCA | 0.9861 / 1.000 / 0.981 | 0.0527 / 0.066 / 0.000 |
| G4 | 4 | Onsager-kT | 0.9936 / 1.000 / 0.947 | 0.0240 / 0.030 / 0.000 |
| G4 | 4 | Onsager-online | 0.9870 / 1.000 / 0.982 | 0.0279 / 0.035 / 0.000 |
| G4 | 8 | SCA | 0.9861 / 1.000 / 0.981 | 0.0572 / 0.072 / 0.000 |
| G4 | 8 | Onsager-kT | 0.9936 / 1.000 / 0.947 | 0.0193 / 0.024 / 0.000 |
| G4 | 8 | Onsager-online | 0.9870 / 1.000 / 0.982 | 0.0286 / 0.036 / 0.000 |
| G5 | 2 | TEC | 0.9887 / 1.000 / 0.891 | 0.0240 / 0.030 / 0.000 |
| G5 | 2 | Onsager-online | 0.9754 / 1.000 / 0.979 | 0.0247 / 0.031 / 0.000 |
| G5 | 4 | TEC | 0.9887 / 1.000 / 0.893 | 0.0371 / 0.047 / 0.000 |
| G5 | 8 | TEC | 0.9887 / 1.000 / 0.893 | 0.0302 / 0.038 / 0.000 |
| G14 | 2 | SCA | 0.9914 / 1.000 / 0.759 | 0.0131 / 0.028 / 0.000 |
| G14 | 4 | SCA | 0.9914 / 1.000 / 0.761 | 0.0323 / 0.069 / 0.000 |
| G14 | 8 | SCA | 0.9914 / 1.000 / 0.760 | 0.0304 / 0.065 / 0.000 |
| G15 | 2 | TEC | 0.9898 / 1.000 / 0.499 | 0.0142 / 0.030 / 0.000 |
| G15 | 2 | Onsager-kT | 0.9880 / 1.000 / 0.557 | 0.0146 / 0.031 / 0.000 |
| G15 | 4 | TEC | 0.9897 / 1.000 / 0.506 | 0.0419 / 0.089 / 0.000 |
| G15 | 4 | Onsager-kT | 0.9880 / 1.000 / 0.562 | 0.0390 / 0.083 / 0.000 |
| G15 | 8 | TEC | 0.9897 / 1.000 / 0.505 | 0.0331 / 0.070 / 0.000 |
| G15 | 8 | Onsager-kT | 0.9881 / 1.000 / 0.564 | 0.0367 / 0.078 / 0.000 |
| G16 | 2 | SCA | 0.9910 / 1.000 / 0.666 | 0.0144 / 0.031 / 0.000 |
| G16 | 4 | SCA | 0.9910 / 1.000 / 0.691 | 0.0352 / 0.076 / 0.000 |
| G16 | 8 | SCA | 0.9910 / 1.000 / 0.690 | 0.0370 / 0.080 / 0.000 |
| G17 | 2 | TEC | 0.9907 / 1.000 / 0.691 | 0.0129 / 0.028 / 0.000 |
| G17 | 2 | Onsager-kT | 0.9909 / 1.000 / 0.769 | 0.0132 / 0.029 / 0.000 |
| G17 | 4 | TEC | 0.9907 / 1.000 / 0.689 | 0.0388 / 0.085 / 0.000 |
| G17 | 4 | Onsager-kT | 0.9910 / 1.000 / 0.767 | 0.0459 / 0.101 / 0.000 |
| G17 | 8 | TEC | 0.9907 / 1.000 / 0.690 | 0.0326 / 0.071 / 0.000 |
| G17 | 8 | Onsager-kT | 0.9909 / 1.000 / 0.767 | 0.0410 / 0.090 / 0.000 |
