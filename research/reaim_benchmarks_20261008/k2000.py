"""K2000 (WK2000_1, the CIM / STATICA / SB / TEC benchmark; best cut 33,337) as a problems.Problem (Max-Cut convention,
J = -w). Source file: fpga/v80_snowball/data/WK2000_1.rud (SHA-256 checked, as research/statica_reproduction_20261003)."""
import hashlib
import math
from pathlib import Path

import numpy as np

import problems as PB

RUD = Path(__file__).resolve().parents[2] / 'fpga' / 'v80_snowball' / 'data' / 'WK2000_1.rud'
RUD_SHA256 = '9ed615e5e18726914f12740b7f9bedb6b69676477ed5958fac42c8ab0c252ba7'
BEST = 33337


def problem():
    raw = RUD.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == RUD_SHA256, 'WK2000_1.rud hash mismatch'
    lines = raw.split(b'\n', 1)
    n, m = map(int, lines[0].split())
    e = np.array(lines[1].split(), dtype=np.int64).reshape(-1, 3)
    assert len(e) == m
    i, j, w = e[:, 0] - 1, e[:, 1] - 1, e[:, 2]
    rows = np.concatenate([i, j]); cols = np.concatenate([j, i]); vals = -np.concatenate([w, w]).astype(np.float64)
    meta = dict(ei=i, ej=j, w=w, W=int(w.sum()), ref=BEST, target=33000, m=len(i))
    return PB.Problem('K2000', 'mcp', n, rows, cols, vals, 0.0, np.zeros(n), math.sqrt(2.0 * len(i) / n), meta)
