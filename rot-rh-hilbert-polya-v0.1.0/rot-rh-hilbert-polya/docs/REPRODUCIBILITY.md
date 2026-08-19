# Reproducibility protocol

## Firewall discipline

A valid reproduction must preserve this ordering:

1. choose all construction parameters;
2. run arithmetic/Gamma/ROT operator construction;
3. write frozen NPZ/JSON artifacts and SHA-256 digest;
4. do not change any construction choice;
5. only then run Xi/zero scoring.

If a parameter is changed because a post-freeze Xi score looked bad, that is a new experiment and must receive a new frozen manifest.

## Precision discipline

For any cancellation-sensitive result:

- repeat with increased quadrature/working precision;
- compare at least two independent numerical formulations where feasible;
- report the large components as well as the small residual;
- test mesh refinement and endpoint stability;
- use arithmetic surrogates/permutations to check whether the observed cancellation is structurally special.

## Heavy-run note

The default release commands are computational research jobs, not unit tests. GitHub Actions checks code integrity and small mathematical invariants only.
