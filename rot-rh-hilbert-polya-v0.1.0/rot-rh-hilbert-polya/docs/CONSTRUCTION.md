# Construction of the current ROT–RH operator

## 1. Arithmetic input

The baseline uses the von Mangoldt function \(\Lambda(n)\), so both primes and prime powers enter naturally. With exponential regulator \(L\), the frozen phase density uses weights of the form

\[
w_n=\frac{\Lambda(n)}{\log n\sqrt n}\,e^{-n/L}.
\]

The phase combines the Archimedean Riemann–Siegel/Gamma contribution with the arithmetic oscillatory contribution. Ordered secular crossings are found **without consulting known zero ordinates**.

## 2. Jacobi realization

For each requested cutoff and depth, the ordered support is converted into a real symmetric Jacobi section. Structural diagnostics include

- self-adjointness defect;
- positivity of off-diagonal Jacobi coefficients;
- eigenvalue preservation;
- orthogonality defect;
- Carleman-type partial sums;
- inverse spectral traces.

The point of the Jacobi realization is to replace a bare list of spectral crossings with an actual finite self-adjoint operator section.

## 3. ROT heat metric

The v97 stage applies the fixed ROT heat parameter `tau` to define the metric used in the Jacobi section. It is part of the **pre-freeze** construction and is not selected by Xi scoring in that run.

## 4. Freeze boundary

The final finite operator arrays are written to a compressed NPZ archive, accompanied by a JSON manifest and SHA-256 digest. After this point the construction is immutable for the scoring pass.

## 5. Infinite-tail completion

A finite section misses high-energy spectral mass. v98–v100.1 model that missing contribution through inverse moments.

For moment index \(m\), with \(p=2m\), the smooth Weyl contribution is

\[
T_m^{\rm Weyl}=
\frac{E_c^{1-2m}}{2\pi(2m-1)}
\left[
\log\frac{E_c}{2\pi}+\frac{1}{2m-1}
\right].
\]

v100.1 adds an exact-minus-leading Archimedean density correction, a frozen Mangoldt oscillatory correction, and the half-level midpoint term

\[
T_m^{\rm midpoint}=-\frac12 E_c^{-2m}.
\]

No Xi values are used to compute these tail coefficients.

## 6. Completed determinant

With positive finite roots \(E_j\),

\[
D_N(t)=\prod_{j=1}^N\left(1-\frac{t^2}{E_j^2}\right).
\]

The moment-completed determinant is

\[
D_M(t)=D_N(t)
\exp\left[-\sum_{m=1}^{M}\frac{T_m t^{2m}}{m}\right].
\]

Only after this function has been frozen is it compared with normalized Xi.
