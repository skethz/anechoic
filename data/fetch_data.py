#!/usr/bin/env python3
"""Download the benchmark instances from their original sources, check their SHA-256 values against the records of the
experiments, and place them where the scripts expect them.

Usage, from the repository root:  python3 data/fetch_data.py [--only k2000|gset|tsplib]"""
import hashlib
import json
import struct
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
K2000_URL = ('https://raw.githubusercontent.com/hariby/SA-complete-graph/'
             '785d664eddb65a5c90a03f5c2253a7cb2cf6ef8b/WK2000_1.rud')
K2000_SHA256 = '9ed615e5e18726914f12740b7f9bedb6b69676477ed5958fac42c8ab0c252ba7'
K2000_BIN_SHA256 = 'a3058db6cffdbb632bc46223958cc7f8dce17da22eb2e4bbed6f60f811e6c0e9'


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, expected, dest):
    dest = Path(dest)
    if dest.exists() and sha256(dest.read_bytes()) == expected:
        print('present ', dest.relative_to(ROOT))
        return dest.read_bytes()
    data = urllib.request.urlopen(url, timeout=120).read()
    if sha256(data) != expected:
        raise SystemExit(f'SHA-256 mismatch for {url}')
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    print('fetched ', dest.relative_to(ROOT))
    return data


def k2000_bin(rud):
    """Bit-packed couplings read by the FPGA host programs and the golden models (sca_ref.hpp, load_graph):
    magic 'SNOWGPU1', uint32 N, then N rows of ceil(N/64) little-endian 64-bit words; bit y of row x is 1 where
    J_xy = +1, that is, where the edge weight is -1 (J = -w). The diagonal bit is 0."""
    lines = rud.decode().split('\n')
    n = int(lines[0].split()[0])
    words = (n + 63) // 64
    rows = [[0] * words for _ in range(n)]
    for line in lines[1:]:
        if not line.strip():
            continue
        i, j, w = map(int, line.split())
        i, j = i - 1, j - 1
        if w == -1:
            rows[i][j >> 6] |= 1 << (j & 63)
            rows[j][i >> 6] |= 1 << (i & 63)
    return b'SNOWGPU1' + struct.pack('<I', n) + b''.join(struct.pack('<%dQ' % words, *r) for r in rows)


def k2000():
    d = ROOT / 'fpga' / 'v80_snowball' / 'data'
    packed = k2000_bin(fetch(K2000_URL, K2000_SHA256, d / 'WK2000_1.rud'))
    if sha256(packed) != K2000_BIN_SHA256:
        raise SystemExit('K2000.bin: conversion does not reproduce the file used in the experiments')
    (d / 'K2000.bin').write_bytes(packed)
    print('built   ', (d / 'K2000.bin').relative_to(ROOT))


def gset():
    base = ROOT / 'research' / 'gset_20261007'
    for entry in json.loads((base / 'manifest.json').read_text())['instances'].values():
        fetch(entry['url'], entry['sha256'], base / entry['file'])


def tsplib():
    base = ROOT / 'research' / 'reaim_benchmarks_20261008'
    for entry in json.loads((base / 'data' / 'tsplib_manifest.json').read_text())['instances'].values():
        fetch(entry['url'], entry['sha256'], base / entry['file'])
        if 'opt_tour_url' in entry:
            fetch(entry['opt_tour_url'], entry['opt_tour_sha256'], base / entry['opt_tour_file'])


if __name__ == '__main__':
    only = sys.argv[2] if len(sys.argv) > 2 and sys.argv[1] == '--only' else None
    for name, step in (('k2000', k2000), ('gset', gset), ('tsplib', tsplib)):
        if only in (None, name):
            step()
