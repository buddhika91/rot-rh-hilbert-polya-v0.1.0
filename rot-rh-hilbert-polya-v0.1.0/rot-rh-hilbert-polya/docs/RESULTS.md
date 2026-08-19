# Current numerical results

## v98 tail-completion benchmark

A reported post-freeze run at final `L=16000`, depth 48 showed

- raw finite determinant symmetric RMS: `1.564988e-02`;
- Weyl-tail completed determinant symmetric RMS: `5.793763e-05`.

This is an improvement by roughly 270×. It is evidence that a substantial part of the finite determinant discrepancy behaves like missing high-energy spectral mass.

## v99 first inverse moment

At depth 48 the reported cutoff was

`Ecut = 139.732482236`.

Representative first-moment terms were

- prefix S1: `1.8452179660e-02`;
- Weyl tail: `+4.6719998137e-03`;
- Archimedean correction: `-8.1021046510e-10`;
- prime-density correction: `+6.3981900359e-06`.

## v100.1 full arithmetic tail

Reported zero-free tail construction:

| depth | Ecut | Weyl | Archimedean | prime | midpoint | completed S1 |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 105.446614 | +5.766178928e-03 | -1.885e-09 | -6.086669216e-06 | -4.496811980e-05 | 2.310493074494e-02 |
| 48 | 139.732482 | +4.671999814e-03 | -8.102e-10 | +6.398190036e-06 | -2.560797613e-05 | 2.310496887752e-02 |

Depth 32→48 convergence:

- moment relative L2: `1.650e-06`;
- determinant relative L2: `6.442e-07`.

These are **depth-convergence metrics of the frozen completed construction**, not a theorem or direct proof-error against Xi.

## Status

`RH_proof = false` remains the correct status. The exact Xi determinant identity and the global cutoff limit are open.
