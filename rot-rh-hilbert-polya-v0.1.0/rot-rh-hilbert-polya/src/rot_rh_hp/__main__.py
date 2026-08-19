from __future__ import annotations
import argparse, json
from pathlib import Path
from .provenance import verify_sources


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser(prog='python -m rot_rh_hp')
    ap.add_argument('command', nargs='?', default='status', choices=['status','verify'])
    args = ap.parse_args()
    root = repo_root()
    if args.command == 'verify':
        results = verify_sources(root)
        print(json.dumps(results, indent=2))
        raise SystemExit(0 if all(results.values()) else 1)
    data = json.loads((root/'results/reported_metrics.json').read_text(encoding='utf-8'))
    print('ROT-RH Hilbert-Polya research release')
    print('Status: numerical candidate; RH proof = FALSE')
    print(f"v98 raw determinant RMS       : {data['v98']['raw_det_symmetric_rms']:.9e}")
    print(f"v98 tail-completed RMS        : {data['v98']['weyl_completed_det_symmetric_rms']:.9e}")
    print(f"v100.1 depth determinant L2   : {data['v100_1']['depth_32_to_48_det_relative_l2']:.9e}")
    print('Exact Fredholm identity       : OPEN')
    print('Global cutoff operator limit  : OPEN')

if __name__ == '__main__':
    main()
