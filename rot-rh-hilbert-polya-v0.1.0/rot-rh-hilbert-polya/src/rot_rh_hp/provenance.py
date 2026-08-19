from __future__ import annotations
import hashlib, json
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify_sources(repo_root: str | Path) -> dict[str, bool]:
    root = Path(repo_root)
    records = json.loads((root/'results/source_provenance.json').read_text(encoding='utf-8'))['files']
    return {r['path']: sha256_file(root/r['path']) == r['sha256'] for r in records}
