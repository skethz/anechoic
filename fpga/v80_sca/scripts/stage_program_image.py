#!/usr/bin/env python3
"""Copy a manifest-verified SCA image to the single image path allowed by the scoped sudo rule. Never programs."""
import argparse
import hashlib
import json
import os
import shutil
import tempfile
import time
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--manifest', type=Path, required=True)
args = parser.parse_args()
root = Path('/scratch/USER/sca_v80_20261003')
manifest_path = args.manifest.resolve(strict=True)
if not manifest_path.is_relative_to(root):
    raise RuntimeError('Manifest must be inside the SCA scratch project')
manifest = json.loads(manifest_path.read_text())
entry = manifest['files']['pdi']
source = Path(entry['path']).resolve(strict=True)
if not source.is_relative_to(root):
    raise RuntimeError('Image must be inside the SCA scratch project')
data = source.read_bytes()
if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
    raise RuntimeError('Image differs from its saved manifest')
destination = Path('/scratch/USER/snowball_v80_20260928/programming')  # fixed by the sudo rule
if destination.is_symlink() or destination.resolve() != destination or not destination.is_dir():
    raise RuntimeError('Programming directory must exist and not redirect elsewhere')
# Preserve the record of what was staged before (its image stays in its own build directory).
previous = destination / 'staged_image.json'
if previous.exists():
    keep = root / 'results' / f'previous_staged_image_{time.strftime("%Y%m%dT%H%M%S")}.json'
    shutil.copy2(previous, keep)
fd, temporary = tempfile.mkstemp(prefix='.image-', dir=destination)
try:
    with os.fdopen(fd, 'wb') as file:
        file.write(data)
        file.flush()
        os.fsync(file.fileno())
    os.replace(temporary, destination / 'snowball_v80.pdi')
finally:
    Path(temporary).unlink(missing_ok=True)
staged = (destination / 'snowball_v80.pdi').read_bytes()
assert hashlib.sha256(staged).hexdigest() == entry['sha256']
previous.write_text(json.dumps({
    'source_manifest': str(manifest_path), 'logic_uuid': manifest['logic_uuid'],
    'sha256': entry['sha256'], 'bytes': len(data), 'project': 'sca_v80_20261003',
    'programmed': False}, indent=2) + '\n')
print(json.dumps({'staged_image': str(destination / 'snowball_v80.pdi'), 'logic_uuid': manifest['logic_uuid'],
                  'sha256': entry['sha256'], 'programmed': False}))
