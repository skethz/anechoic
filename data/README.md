# Data

The benchmark instances are not redistributed here. From the repository root, `python3 data/fetch_data.py` downloads them from their original sources. It checks every file against the SHA-256 recorded when the experiments were run, and places it where the scripts expect it.

| Data | Source | Placed in | Record |
|---|---|---|---|
| K2000 (`WK2000_1.rud`): complete graph, 2,000 nodes, ±1 weights | [hariby/SA-complete-graph](https://github.com/hariby/SA-complete-graph), branch `WK2000`, commit `785d664` (MIT license) | `fpga/v80_snowball/data/` | SHA-256 in `fetch_data.py` |
| `K2000.bin`: bit-packed couplings for the FPGA host programs and the golden models | Converted from `WK2000_1.rud` by `fetch_data.py` | `fpga/v80_snowball/data/` | The conversion is checked against the SHA-256 of the file used in the experiments |
| G-set: 51 instances (G1–G47, G51–G54) | [Y. Ye, Stanford University](https://web.stanford.edu/~yyye/yyye/Gset/) | `research/gset_20261007/data/` | `research/gset_20261007/manifest.json` |
| TSPLIB95: gr17, gr21, gr24, fri26, bayg29, bays29, and the optimal tours | [TSPLIB95, Heidelberg University](https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) | `research/reaim_benchmarks_20261008/data/tsplib/` | `research/reaim_benchmarks_20261008/data/tsplib_manifest.json` |

On the experiment hosts, the FPGA and GPU run scripts read the instances from a per-task data directory, for example `$TASK_ROOT/data/K2000.bin` or `/scratch/USER/.../data/gset/G22`. Copy the downloaded files there.

The authors computed the graph-partitioning reference partitions (`research/reaim_benchmarks_20261008/data/gpp_reference.json`) by simulated annealing with balance-preserving swaps. They are included.
