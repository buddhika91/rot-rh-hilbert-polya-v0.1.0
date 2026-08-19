# ROT–RH Hilbert–Pólya Operator

> **Research status:** numerical/operator-construction program; **not a proof of the Riemann Hypothesis**.

This repository packages the current ROT–RH Hilbert–Pólya candidate centered on the **v97 → v100.1 frozen arithmetic/Fredholm pipeline**, with the later v102.1 finite-part extension and v129 red-team audit retained as research branches.

The central experiment is deliberately split into two phases:

1. **Pre-freeze construction:** build the spectral operator and arithmetic tail from Gamma/Archimedean structure, Mangoldt-weighted prime-power data, the fixed ROT heat metric, and cutoff/depth rules. No known zeta-zero ordinates and no `Xi(t)` values are used to select or tune the construction.
2. **Post-freeze scoring:** hash/freeze the constructed object, then evaluate the normalized Riemann function and compare the frozen spectral determinant with `Xi(t)/Xi(0)`.

The exact theorem target is

\[
D_H(t)=\det_{\mathrm{reg}}\!\left(I-t^2H^{-2}\right)=\frac{\Xi(t)}{\Xi(0)},
\]

for a cutoff-independent infinite self-adjoint operator \(H\). The current code provides strong finite-cutoff numerical evidence and a fixed-\(L\) operator construction; the global cutoff limit and exact entire-function identity remain open.

## Current reported benchmark

The release records the latest reported v97–v100.1 benchmark without pretending to have recomputed the heavy jobs inside CI:

| Quantity | Reported value |
|---|---:|
| v98 raw finite-determinant symmetric RMS | `1.564988e-02` |
| v98 Weyl-tail completed symmetric RMS | `5.793763e-05` |
| Improvement factor | about `270×` |
| v100.1 depth 32→48 moment relative L2 | `1.650e-06` |
| v100.1 depth 32→48 determinant relative L2 | `6.442e-07` |
| v100.1 completed S1, depth 32 | `2.310493074494e-02` |
| v100.1 completed S1, depth 48 | `2.310496887752e-02` |
| v100.1 pre-freeze Xi/zeta/known-zero access | `NONE` |
| RH proof | `FALSE` |

The depth-convergence number is an **internal convergence diagnostic**, not an error bar proving equality to \(\Xi\).

## Repository map

```text
rot-rh-hilbert-polya/
├── research/
│   ├── legacy_pipeline/          # Original recovered v97–v102.1 research scripts
│   └── red_team/                 # v129 adaptive-mesh surrogate audit
├── src/rot_rh_hp/                # Small reusable reference library
├── tests/                        # Fast integrity, firewall, determinant and tail tests
├── docs/                         # Theory, construction, results, limitations, roadmap
├── results/                      # Machine-readable reported benchmark + provenance
├── manuscript/                   # Article and paper outline
├── scripts/                      # PowerShell/bash reproduction helpers
└── .github/workflows/ci.yml      # GitHub Actions smoke CI
```

## Installation

Python 3.10+ is supported; Python 3.11 is recommended.

### Windows PowerShell

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -e ".[dev]"
pytest
python -m rot_rh_hp status
```

### Linux/macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e '.[dev]'
pytest
python -m rot_rh_hp status
```

## Reproduce the release pipeline

The original research scripts are kept with their historical filenames because later stages dynamically import earlier ones.

### Windows

```powershell
.\scripts\run_v97_v100.ps1
```

### Linux/macOS

```bash
bash scripts/run_v97_v100.sh
```

These are intentionally **heavy numerical runs**. CI only performs smoke/integrity tests.

## Exact construction chain

At finite exponential cutoff \(L\), v97 builds ordered arithmetic roots and a symmetric Jacobi realization. The construction is based on a Gamma/Archimedean phase plus exponentially regulated Mangoldt data and a ROT heat metric. It then freezes the finite operator data before any zeta/zero scoring.

The v100.1 tail law writes the missing inverse moments schematically as

\[
T_m=T_m^{\mathrm{Weyl}}+T_m^{\mathrm{Arch}}+T_m^{\mathrm{prime}}+T_m^{\mathrm{midpoint}},
\]

and completes the finite determinant by

\[
D(t)=D_N(t)\exp\left[-\sum_{m=1}^{M}\frac{T_m t^{2m}}{m}\right].
\]

See [`docs/CONSTRUCTION.md`](docs/CONSTRUCTION.md) and [`docs/FREDHOLM_IDENTITY.md`](docs/FREDHOLM_IDENTITY.md).

## Scientific claim boundary

This repository supports the following conservative statement:

> A zero-free arithmetic/ROT construction produces self-adjoint finite/fixed-cutoff spectral sections and, after an independently derived arithmetic tail completion, exhibits strong post-freeze agreement and depth convergence against the normalized Riemann Xi determinant.

It **does not** establish:

- convergence of the full \(L\to\infty\) operator sequence;
- a complete domain/resolvent theorem for the global limiting operator;
- uniform compact-set convergence of the regularized determinants;
- the exact entire-function identity with \(\Xi\);
- the Riemann Hypothesis.

Those are the theorem-level bottlenecks documented in [`docs/ROADMAP.md`](docs/ROADMAP.md).

## Provenance

The historical Python scripts were recovered from the project library and copied byte-for-byte into `research/`. Their SHA-256 hashes are recorded in `results/source_provenance.json` and checked by the test suite.

## License

MIT for the repository scaffolding and reusable reference code. Historical research scripts are included as project source under the same repository license unless superseded by a future explicit notice.
