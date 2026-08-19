# The Current ROT–RH Hilbert–Pólya Operator: Construction, Fredholm Closure, and Numerical Performance

## Abstract

The ROT–RH program seeks a Hilbert–Pólya realization of the Riemann spectrum from arithmetic and recursive structure rather than from a supplied list of zeta zeros. Its current v97–v100.1 architecture constructs a self-adjoint finite/fixed-cutoff Jacobi operator from Gamma/Archimedean structure, Mangoldt-weighted prime powers and a fixed ROT heat metric, freezes the resulting operator, derives a zero-free arithmetic completion of the unresolved spectral tail, and only then compares the frozen Fredholm determinant with `Xi(t)/Xi(0)`.

The key result is methodological and spectral rather than merely a low zero RMSE. A reported v98 tail completion lowered a representative symmetric determinant RMS from about `1.56e-2` to `5.79e-5`. The v100.1 full arithmetic completion then showed depth 32→48 internal determinant convergence at about `6.44e-7` and moment convergence at about `1.65e-6`. The tail decomposes into Weyl, Archimedean, frozen-prime and midpoint terms.

These results make the construction a nontrivial numerical Hilbert–Pólya candidate, but not an RH proof. The remaining theorem-level problem is to prove that the cutoff family converges to an infinite self-adjoint operator and that its regularized determinant is exactly the normalized Riemann Xi function as an entire-function identity.

## Construction

The arithmetic phase uses Mangoldt weights

\[
w_n=\Lambda(n)(\log n\sqrt n)^{-1}e^{-n/L},
\]

combined with the Archimedean phase. Ordered secular roots are converted into a symmetric Jacobi realization. ROT enters through the fixed heat-metric construction used in the v97 section.

## Fredholm completion

For finite positive roots \(E_j\),

\[
D_N(t)=\prod_j(1-t^2/E_j^2).
\]

The missing tail inverse moments are modeled by

\[
T_m=T_m^{W}+T_m^{A}+T_m^{P}+T_m^{M},
\]

and the determinant is completed as

\[
D(t)=D_N(t)\exp[-\sum_m T_m t^{2m}/m].
\]

The complete pre-freeze object is hashed before Xi is evaluated.

## Exact target

The final desired theorem is

\[
\det_{reg}(I-t^2H^{-2})=\Xi(t)/\Xi(0).
\]

A proof of this identity for a self-adjoint limiting \(H\) would be far stronger than matching any finite number of zeros. It would identify the entire spectral determinant with Xi.

## Present status

The fixed-cutoff program is computationally mature enough to support serious convergence and falsification tests. The global `L→∞` operator limit, determinant-class theorem, uniform compact-set convergence and exact Xi identity remain open. Those are the correct next targets.
