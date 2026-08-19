# Relevance to the exact Fredholm identity

The numerical target is not merely a root fit. The intended theorem is an operator identity.

Let \(H\) denote the hoped-for infinite self-adjoint Hilbert–Pólya operator with symmetric spectrum \(\{\pm E_n\}\). Formally,

\[
\det\!\left(I-t^2H^{-2}\right)
=\prod_n\left(1-\frac{t^2}{E_n^2}\right).
\]

The ROT–RH target is

\[
\boxed{\det_{\rm reg}(I-t^2H^{-2})=\Xi(t)/\Xi(0).}
\]

Why this matters:

1. a self-adjoint \(H\) has real spectrum;
2. zeros of its determinant occur at spectral values;
3. equality with normalized Xi would identify the entire zero set, including multiplicity;
4. with all analytic details proved, this would provide the Hilbert–Pólya bridge to RH.

## What the numerics test

The current finite pipeline tests whether

\[
D_{L,N,M}(t)\to \Xi(t)/\Xi(0)
\]

as depth, arithmetic cutoff and tail order are increased.

The key methodological improvement is that `D_{L,N,M}` is frozen *before* Xi is evaluated.

## What is still missing

Numerical agreement, regardless of precision, is not the exact identity. A proof requires at least:

- a well-defined cutoff-independent limiting operator;
- self-adjointness of that limit on an explicit dense domain;
- a justified determinant class/regularization;
- compact-set convergence of finite/regularized determinants;
- analytic identification of the limiting entire function with Xi.

The present repository keeps these items explicitly open.
