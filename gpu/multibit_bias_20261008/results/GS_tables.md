## Amendment 5: GPU-selected G-set schedules

Verification: 1996/1996 traced trials bit-exact against run_trial_bias in 499 runs (one per unique selected cohort); device/host mismatches 0; final-field mismatches 0.

Held-out cohorts: 499/499 measured, 128 batches each; device/host cut mismatches 0.


### Original grid: GPU (own selection) against V80 (own selection, A3 board)

| Rule | GPU primary (ms) | V80 primary (ms) | V80/GPU | GPU faster (of 51) | GPU secondary (ms) | V80 secondary, A3 at E12 (ms) | V80/GPU | GPU faster | V80 secondary, own selection, modelled (ms) | V80/GPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SCA | 0.5938 | 0.2411 | 0.41 | 6 | 0.08599 | 0.18532 | 2.16 | 47 | 0.09855 | 1.15 |
| TEC | 0.5020 | 0.1883 | 0.38 | 4 | 0.06650 | 0.14316 | 2.15 | 46 | 0.08623 | 1.30 |
| Onsager-kT | 0.3462 | 0.1211 | 0.35 | 5 | 0.04116 | 0.09070 | 2.20 | 50 | 0.05810 | 1.41 |
| Onsager-online | 0.5882 | 0.1907 | 0.32 | 9 | 0.07008 | 0.14669 | 2.09 | 50 | 0.08807 | 1.26 |
| best of 4 rules | 0.3406 | 0.1145 | 0.34 | 5 | 0.04026 | 0.07782 | 1.93 | 50 | 0.05287 | 1.31 |

Per graph class (original grid; geometric means in ms; V80/GPU in brackets; GPU wins of n):

| Rule | Class | n | GPU primary | V80 primary | GPU wins | GPU secondary | V80 secondary A3 | GPU wins | V80 secondary modelled |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SCA | random | 25 | 0.2133 | 0.0755 [0.35] | 0 | 0.01496 | 0.05042 [3.37] | 25 | 0.01886 [1.26] |
| SCA | toroidal | 6 | 5.8804 | 6.0696 [1.03] | 3 | 5.16652 | 6.16648 [1.19] | 4 | 5.71157 [1.11] |
| SCA | planar | 20 | 1.0734 | 0.3908 [0.36] | 3 | 0.22395 | 0.32958 [1.47] | 18 | 0.23030 [1.03] |
| TEC | random | 25 | 0.2075 | 0.0529 [0.26] | 0 | 0.01194 | 0.03519 [2.95] | 25 | 0.01586 [1.33] |
| TEC | toroidal | 6 | 6.4559 | 6.1188 [0.95] | 2 | 6.22535 | 6.15425 [0.99] | 2 | 6.85190 [1.10] |
| TEC | planar | 20 | 0.7037 | 0.3238 [0.46] | 2 | 0.14575 | 0.26762 [1.84] | 19 | 0.19263 [1.32] |
| Onsager-kT | random | 25 | 0.1845 | 0.0374 [0.20] | 0 | 0.00889 | 0.02380 [2.68] | 25 | 0.01171 [1.32] |
| Onsager-kT | toroidal | 6 | 1.8315 | 1.4093 [0.77] | 1 | 1.48456 | 1.42282 [0.96] | 5 | 1.43934 [0.97] |
| Onsager-kT | planar | 20 | 0.4611 | 0.2522 [0.55] | 4 | 0.09536 | 0.21153 [2.22] | 20 | 0.16420 [1.72] |
| Onsager-online | random | 25 | 0.2895 | 0.0506 [0.17] | 0 | 0.01349 | 0.03480 [2.58] | 25 | 0.01506 [1.12] |
| Onsager-online | toroidal | 6 | 4.0740 | 4.3112 [1.06] | 4 | 3.63192 | 4.35897 [1.20] | 5 | 5.03025 [1.39] |
| Onsager-online | planar | 20 | 0.7983 | 0.3926 [0.49] | 5 | 0.16822 | 0.32024 [1.90] | 20 | 0.23793 [1.41] |
| best of 4 rules | random | 25 | 0.1834 | 0.0354 [0.19] | 0 | 0.00856 | 0.01977 [2.31] | 25 | 0.01039 [1.21] |
| best of 4 rules | toroidal | 6 | 1.8315 | 1.4093 [0.77] | 1 | 1.48456 | 1.42282 [0.96] | 5 | 1.43934 [0.97] |
| best of 4 rules | planar | 20 | 0.4456 | 0.2339 [0.52] | 4 | 0.09449 | 0.18041 [1.91] | 20 | 0.14991 [1.59] |

### Extended grid: GPU (own selection) against V80 (own selection, A3 board)

| Rule | GPU primary (ms) | V80 primary (ms) | V80/GPU | GPU faster (of 51) | GPU secondary (ms) | V80 secondary, A3 at E12 (ms) | V80/GPU | GPU faster | V80 secondary, own selection, modelled (ms) | V80/GPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SCA | 0.5970 | 0.2396 | 0.40 | 6 | 0.08555 | 0.17881 | 2.09 | 49 | 0.10226 | 1.20 |
| TEC | 0.4600 | 0.1856 | 0.40 | 7 | 0.06190 | 0.13040 | 2.11 | 49 | 0.08139 | 1.31 |
| Onsager-kT | 0.3512 | 0.1106 | 0.31 | 4 | 0.03821 | 0.07398 | 1.94 | 49 | 0.04742 | 1.24 |
| Onsager-online | 0.6247 | 0.1897 | 0.30 | 7 | 0.07434 | 0.14544 | 1.96 | 48 | 0.08863 | 1.19 |
| best of 4 rules | 0.3445 | 0.1043 | 0.30 | 5 | 0.03740 | 0.06583 | 1.76 | 49 | 0.04453 | 1.19 |

Per graph class (extended grid; geometric means in ms; V80/GPU in brackets; GPU wins of n):

| Rule | Class | n | GPU primary | V80 primary | GPU wins | GPU secondary | V80 secondary A3 | GPU wins | V80 secondary modelled |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SCA | random | 25 | 0.2133 | 0.0730 [0.34] | 0 | 0.01500 | 0.04645 [3.10] | 25 | 0.01885 [1.26] |
| SCA | toroidal | 6 | 5.8260 | 7.4092 [1.27] | 4 | 5.03944 | 7.29311 [1.45] | 6 | 8.58671 [1.70] |
| SCA | planar | 20 | 1.0911 | 0.3782 [0.35] | 2 | 0.22202 | 0.31695 [1.43] | 18 | 0.22413 [1.01] |
| TEC | random | 25 | 0.2074 | 0.0491 [0.24] | 0 | 0.01117 | 0.02848 [2.55] | 25 | 0.01413 [1.26] |
| TEC | toroidal | 6 | 5.4843 | 6.2156 [1.13] | 3 | 5.19030 | 6.20887 [1.20] | 4 | 6.32187 [1.22] |
| TEC | planar | 20 | 0.5919 | 0.3410 [0.58] | 4 | 0.13937 | 0.27410 [1.97] | 20 | 0.19688 [1.41] |
| Onsager-kT | random | 25 | 0.1845 | 0.0389 [0.21] | 0 | 0.00817 | 0.02006 [2.46] | 25 | 0.01131 [1.38] |
| Onsager-kT | toroidal | 6 | 1.9466 | 1.6238 [0.83] | 1 | 1.63550 | 1.61377 [0.99] | 4 | 1.77866 [1.09] |
| Onsager-kT | planar | 20 | 0.4696 | 0.1820 [0.39] | 3 | 0.08515 | 0.14998 [1.76] | 20 | 0.09596 [1.13] |
| Onsager-online | random | 25 | 0.2895 | 0.0506 [0.17] | 0 | 0.01349 | 0.03480 [2.58] | 25 | 0.01542 [1.14] |
| Onsager-online | toroidal | 6 | 4.9147 | 5.1142 [1.04] | 3 | 4.34877 | 5.21076 [1.20] | 4 | 6.19914 [1.43] |
| Onsager-online | planar | 20 | 0.8800 | 0.3684 [0.42] | 4 | 0.18520 | 0.29701 [1.60] | 19 | 0.22064 [1.19] |
| best of 4 rules | random | 25 | 0.1834 | 0.0362 [0.20] | 0 | 0.00787 | 0.01667 [2.12] | 25 | 0.01021 [1.30] |
| best of 4 rules | toroidal | 6 | 1.9466 | 1.6238 [0.83] | 1 | 1.63550 | 1.61377 [0.99] | 4 | 1.77866 [1.09] |
| best of 4 rules | planar | 20 | 0.4505 | 0.1722 [0.38] | 4 | 0.08451 | 0.14039 [1.66] | 20 | 0.09279 [1.10] |

### For reference: the earlier identical-schedule comparison (GPU on the study's extended-grid E12 configurations of phase G, against the V80's extended-grid A3 cohorts, which are the same configurations)

| Rule | GPU primary, best B per instance (ms) | GPU primary, B = 60 (ms) | V80 primary (ms) | V80/GPU (best B; B = 60) | GPU secondary, best B (ms) | GPU secondary, B = 240 (ms) | V80 secondary (ms) | V80/GPU (best B; B = 240) |
|---|---:|---:|---:|---|---:|---:|---:|---|
| SCA | 0.8154 | 0.8316 | 0.2396 | 0.29; 0.29 | 0.13622 | 0.13707 | 0.17881 | 1.31; 1.30 |
| TEC | 0.5846 | 0.5938 | 0.1856 | 0.32; 0.31 | 0.09035 | 0.09095 | 0.13040 | 1.44; 1.43 |
| Onsager-kT | 0.4403 | 0.4438 | 0.1106 | 0.25; 0.25 | 0.05647 | 0.05672 | 0.07398 | 1.31; 1.30 |
| Onsager-online | 0.8863 | 0.9182 | 0.1897 | 0.21; 0.21 | 0.12143 | 0.13883 | 0.14544 | 1.20; 1.05 |

### GPU selections (original grid): distribution of (S, B)

- primary: S=250 B=60: 95, S=250 B=120: 8, S=500 B=60: 25, S=500 B=120: 10, S=500 B=240: 2, S=1000 B=60: 19, S=1000 B=120: 12, S=1000 B=240: 1, S=2000 B=60: 15, S=2000 B=120: 6, S=2000 B=240: 2, S=4000 B=60: 3, S=4000 B=120: 6
- secondary: S=250 B=120: 1, S=250 B=240: 2, S=500 B=120: 5, S=500 B=240: 26, S=1000 B=120: 24, S=1000 B=240: 44, S=2000 B=120: 12, S=2000 B=240: 54, S=4000 B=120: 27, S=4000 B=240: 9
- held-out / predicted TTS, median over labels: primary 0.999, secondary 1.043

### GPU board power at the new schedules (nvidia-smi on the UUID, 200 ms; 40 s per combination)

| Combination (B, variant, class) | samples | board W | module W |
|---|---:|---:|---:|
| B60_R2G1_const | 153 | 301.1 | 398.4 |
| B240_R1G2_const | 153 | 365.1 | 466.0 |
| B60_R2G1_ons | 152 | 270.4 | 371.3 |
| B120_R1G1_ons | 152 | 354.7 | 455.6 |
| B120_R1G1_const | 153 | 361.8 | 499.4 |
| B120_R2G2_const | 153 | 355.9 | 606.5 |
| B240_R2G4_ons | 153 | 354.1 | 606.7 |
| B120_R2G2_ons | 154 | 354.7 | 580.9 |
| idle | 133 | 99.1 | 192.9 |

Energy to solution, geometric means over the 51 instances (GPU: board power of the cohort's combination × TTS; V80: mean A3 board power of its G-set loads, 93.9 W, × TTS):

| Grid | Rule | GPU primary (mJ) | V80 primary (mJ) | GPU secondary (mJ) | V80 secondary A3 (mJ) |
|---|---|---:|---:|---:|---:|
| original | SCA | 186.5 | 22.65 | 31.13 | 17.409 |
| original | TEC | 159.0 | 17.69 | 24.19 | 13.449 |
| original | Onsager-kT | 106.1 | 11.38 | 15.00 | 8.521 |
| original | Onsager-online | 173.2 | 17.91 | 24.85 | 13.780 |
| extended | SCA | 187.5 | 22.51 | 30.92 | 16.798 |
| extended | TEC | 144.5 | 17.44 | 22.52 | 12.250 |
| extended | Onsager-kT | 108.0 | 10.39 | 13.91 | 6.950 |
| extended | Onsager-online | 184.9 | 17.82 | 26.36 | 13.663 |
