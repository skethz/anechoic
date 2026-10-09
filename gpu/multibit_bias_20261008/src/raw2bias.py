#!/usr/bin/env python3
"""Convert a raw little-endian int32 bias array (N values, as written by research/reaim_benchmarks_20261008/problems.py
write_bias) into the SCABIAS1 file read by sca_gpu_bias / check_ref_bias (8-byte magic, uint32 N, uint32 0, N int32).
N is taken from the SCAJINT8 header of the matching coupling file. Values are copied unchanged.
Usage: raw2bias.py J.jint8 bias.raw out.bias"""
import struct, sys
jf, bf, of = sys.argv[1:4]
with open(jf, 'rb') as f:
    hdr = f.read(16)
assert hdr[:8] == b'SCAJINT8', 'not a SCAJINT8 file'
n = struct.unpack('<I', hdr[8:12])[0]
raw = open(bf, 'rb').read()
assert len(raw) == 4 * n, f'bias file has {len(raw)} bytes, expected {4 * n}'
vals = struct.unpack(f'<{n}i', raw)
assert all(-32767 <= v <= 32767 for v in vals), 'bias outside |b| <= 32767'
with open(of, 'wb') as f:
    f.write(b'SCABIAS1'); f.write(struct.pack('<II', n, 0)); f.write(raw)
print(n, 'biases, max |b| =', max(abs(v) for v in vals))
