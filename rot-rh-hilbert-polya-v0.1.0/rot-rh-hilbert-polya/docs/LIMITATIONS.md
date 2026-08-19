# Limitations and anti-overclaim rules

1. **No finite RMSE establishes RH.** A finite-grid comparison can falsify a candidate but cannot prove an entire-function identity.
2. **Self-adjoint finite matrices are not enough.** The essential object is a well-defined infinite cutoff-independent operator.
3. **Tail asymptotics need rigorous error bounds.** v100.1 uses stable asymptotic/integral evaluations, but the theorem target needs uniform remainder control.
4. **The Xi firewall is methodological, not itself a proof.** It protects against direct target fitting, but arithmetic formulae may still encode classical zeta structure.
5. **Red-team failures count.** Earlier interpolation/quadrature artifacts motivate independent numerical methods, precision doubling, surrogate controls and mesh checks.
6. **v129 is not self-contained.** It expects a v123 `_CELLS.csv` artifact not recovered into this release. It is retained for provenance and future continuation, not advertised as a one-command reproducible stage.
