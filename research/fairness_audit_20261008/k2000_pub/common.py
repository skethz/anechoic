"""Path shim for running the frozen algorithm-comparison modules outside the laptop (no original file is edited).
The frozen modules insert their sibling folders into sys.path; this shim also appends them, so that imports resolve
when the modules run from another folder (their SHA-256 values are checked in pub_methods.selftest). AUDIT_ROOT (on gpu-host: /scratch/USER/anechoic_fairness_20261008) mirrors the project's research/ layout."""
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get('AUDIT_ROOT', Path(__file__).resolve().parents[3]))
for sub in ('research/theory_ideas_20261003', 'research/statica_reproduction_20261003', 'research/ablation_20261005',
            'research/reaim_reproduction_20261003', 'research/algorithm_compare_20261007'):
    p = str(ROOT / sub)
    if p not in sys.path:
        sys.path.append(p)
import methods as M  # noqa: E402,F401
