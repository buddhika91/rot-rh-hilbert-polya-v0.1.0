#!/usr/bin/env python3
from pathlib import Path
import hashlib, sys
root=Path(__file__).resolve().parents[1]
failed=[]
for line in (root/'MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
    if not line.strip(): continue
    expected, rel=line.split('  ',1)
    p=root/rel
    if not p.exists():
        failed.append((rel,'missing'))
        continue
    got=hashlib.sha256(p.read_bytes()).hexdigest()
    if got != expected: failed.append((rel,got))
if failed:
    for x in failed: print('FAIL',*x)
    sys.exit(1)
print('Manifest OK')
