# Anechoic

Code and summary data for **"Anechoic: Onsager-Corrected Parallel Annealing for Sub-0.1 ms Dense Max-Cut on an FPGA"** by Seungki Hong (ETH Zurich), Kyeongwon Jeong (Yonsei University) and Taekwang Jang (ETH Zurich). The preprint link will be added here; the LaTeX source is in [`paper/`](paper/).

Stochastic cellular automata (SCA) update all spins of an Ising machine at once. This makes each spin's previous state return to it through its neighbours, an echo that undoes about two-thirds of all flips on K2000. Anechoic subtracts this echo from every decision. Dynamical mean-field theory gives the coefficient as the number of spins with a fractional flip probability, divided by twice the temperature. The engine counts it online (Onsager-online) or scales it with the temperature (Onsager-κT). This repository holds:
- the RTL engine for the AMD Alveo V80;
- the bit-exact GPU kernel;
- the RTL of the ASIC variant;
- the golden models;
- the theory and algorithm studies and the machine-checked proofs;
- the summary data behind the paper's tables and figures.

## Layout

| Folder | Contents |
|---|---|
| `fpga/v80_sca/` | V80 engine. RTL in `src/v6/` (`sca_core*.v`); HLS kernels, host programs and golden models in `src/`; build, run and analysis scripts in `scripts/`. Also the frozen protocols (`PROTOCOL_HW*.md`), the design notes (`DESIGN.md`, `MULTIBIT_SPEC.md`) and the result summaries (`RESULTS_HW.md`, `RESULTS_MULTIBIT.md`, `results/`). |
| `gpu/onsager_v2_20261007/` | Bit-exact kernel for the NVIDIA GH200 (`src/sca_gpu.cu`), with protocol, report and result tables. |
| `gpu/multibit_bias_20261008/` | GPU kernel with multi-bit couplings and bias, with its G-set results. |
| `asic/gf22_engine_v64_1p25_20261007/` | ASIC variant of the v6 engine at 1.25 GHz. RTL in `rtl/`; source RTL and golden model in `orig/`; RTL regression in `sim/` and `scripts/`; derived 12-engine totals in `results/`. |
| `asic/gf22_engine_v64_20261007/` | The same, 1 GHz variant. |
| `asic/gf22_engine_v64_mb_20261008/` | Multi-bit ASIC variant (K = 2 and K = 4) with bias. |
| `research/dmft_sca_20261004/` | Dynamical mean-field theory of synchronous SCA. |
| `research/formal_verification_20261007/` | Lean 4/Mathlib (`SnowballFormal/`) and Rocq (`coq/`) proofs. |
| `research/algorithm_compare_20261007/`, `ablation_20261005/`, `theory_ideas_20261003/`, `iamp_evaluation_20261003/`, `statica_reproduction_20261003/`, `reaim_reproduction_20261003/` | The corrected update rules, the re-implemented baselines and their studies. |
| `research/fairness_audit_20261008/` | Comparison at the baselines' published settings (`k2000_pub/`), G-set speed-ups (`gset/`) and hardware tables (`hw_table2/`). |
| `research/v6_schedule_20261006/` | Schedule selection for the board configurations. |
| `research/gset_20261007/` | G-set study with the engine's model. |
| `research/reaim_benchmarks_20261008/` | Max-Cut, graph partitioning and TSP benchmarks in software. |
| `research/optimum_mitigation_20261007/` | Study at the best-known K2000 cut. |
| `paper/` | LaTeX source of the preprint. |
| `data/` | Download and verification of the benchmark instances. |

## Where the paper's results come from

| Paper | Source |
|---|---|
| Figure 1 (theory and simulation) | `research/dmft_sca_20261004/` (`dmft.py`, `dmft_torch.py`, `plot_theory_vs_sim_v2.py`) |
| Table 1 (machine-checked statements) | `research/formal_verification_20261007/` (`SPEC.md`) |
| Figure 2 (comparison at published settings) | `research/fairness_audit_20261008/k2000_pub/` (`run_pub.py`, `make_fig.py`, `pub_methods.py`), with the algorithm folders above |
| Stability predictions for the correction | `research/theory_ideas_20261003/` (`a2.py`, `PROTOCOL_A2.md`) |
| Figure 3 and the engine design | `fpga/v80_sca/src/v6/`, `fpga/v80_sca/DESIGN.md` |
| Tables 2 and 3 (board TTS99, design evolution) | `fpga/v80_sca/` (`RESULTS_HW.md`, `RESULTS_MULTIBIT.md`, `PROTOCOL_HW*.md`, `results/`) |
| Figure 4 (flips per step) | `research/algorithm_compare_20261007/flips_profile_v2.py` |
| Tables 4 and 5 (Anechoic on FPGA, GPU and ASIC) | `fpga/v80_sca/`, `gpu/onsager_v2_20261007/` (`REPORT.md`), `asic/gf22_engine_v64_1p25_20261007/results/` |
| Table 6 (FPGA resources) | `fpga/v80_sca/results/util_20261007/` |
| Table 7 and the ASIC results | `asic/` (RTL, regression, derived totals). The layout figure is in `paper/`. |
| Table 8 (G-set, measured) | `fpga/v80_sca/RESULTS_MULTIBIT.md` and `PROTOCOL_HW_MB_A3*`, `research/gset_20261007/`, `research/fairness_audit_20261008/gset/` |
| Table 9 (software benchmarks) | `research/reaim_benchmarks_20261008/` (`run_a5*.py`, `analyze_a5.py`, `compact_table*.md`) |
| Comparison with GbSB at the best-known cut | `research/optimum_mitigation_20261007/` |

## Requirements

- **Python 3.11 or newer** with numpy, numba, scipy and matplotlib. `dmft_torch.py` also needs PyTorch. The Neal reference runs need dimod and dwave-samplers. Exact versions are in the protocol files.
- **A C++17 compiler** for the golden models, the vector generators and the host programs.
- **For the GPU kernels:** CUDA 13 and an NVIDIA GH200 (sm_90).
- **For the V80 builds:** AMD Vivado/Vitis 2025.1, the [AVED](https://github.com/Xilinx/AVED) shell and an AMD Alveo V80.
- **For the RTL regressions:** a Verilog simulator (Siemens Questa or AMD Vivado xsim).
- **For the proofs:** Lean 4, at the version in `SnowballFormal/lean-toolchain`, with Mathlib; and the Rocq Prover 8.20.

Run `python3 data/fetch_data.py` once to obtain the benchmark instances; see [`data/README.md`](data/README.md).

## Notes

- **Placeholders.** The experiments ran on the authors' servers. For publication, user names, host names, local paths and device identifiers were replaced by placeholders: `/scratch/USER`, `gpu-host`, `fpga-host` and `GPU-xxxxxxxx-…`. Set them for your environment.
- **Hash records.** The SHA-256 records of the frozen protocols and inputs (`*.sha256`, `SHA256SUMS*.txt`) refer to the unredacted originals. The self-test in `research/fairness_audit_20261008/k2000_pub/pub_methods.py` checks the published files.
- **Raw data.** Per-trial records, spin states and traces are not included because of their size. The summaries, and the scripts that produced them, are included.
- **ASIC.** The ASIC folders contain the RTL, its regression against the golden model and the derived totals. Synthesis and place-and-route scripts and reports are not included, because they depend on a foundry design kit under a non-disclosure agreement.
- **FPGA platform.** The V80 builds use AMD's AVED shell, together with the platform-integration and programming scripts of the authors' earlier V80 design (`fpga/v80_snowball/`). Neither is part of this repository.

## Citation

The preprint entry will be added here.

## License

To be added.
