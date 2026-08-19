# ROT–RH Hilbert–Pólya Operator

## Zero-free arithmetic spectral construction toward an exact Fredholm representation of the Riemann Xi function

> **Research status:** active numerical and operator-theoretic program.  
> **This repository does not claim a proof of the Riemann Hypothesis.**

This repository contains the current **Recursive Observation Theory–Riemann Hypothesis (ROT–RH) Hilbert–Pólya operator program**.

The goal is to construct an infinite self-adjoint operator `H` whose positive spectrum reproduces the ordinates of the non-trivial zeros of the Riemann zeta function,

```math
\rho_n=\frac12+i\gamma_n,
```

through

```math
H\psi_n=\gamma_n\psi_n.
```

The stronger target is the exact Fredholm identity

```math
\boxed{
D_H(t)
=
\det_{\mathrm{reg}}\!\left(I-t^2H^{-2}\right)
=
\frac{\Xi(t)}{\Xi(0)}.
}
```

The present codebase goes beyond direct zero fitting. It contains a zero-free arithmetic construction, a self-adjoint Jacobi realization, a frozen pre-Xi operator object, arithmetic completion of the unresolved spectral tail, post-freeze Xi scoring, convergence diagnostics, and adversarial numerical audits.

The core methodological rule is:

```text
construct first
freeze second
compare with Xi only afterward
```

Known zero ordinates and evaluations of `Xi` are therefore treated as **validation data**, not construction data.

---

# 1. Mathematical target

For

```math
\Re(s)>1,
```

the Riemann zeta function is

```math
\zeta(s)
=
\sum_{n=1}^{\infty}\frac1{n^s}
=
\prod_p\left(1-p^{-s}\right)^{-1}.
```

The Riemann Hypothesis states that all non-trivial zeros satisfy

```math
\boxed{
\rho_n=\frac12+i\gamma_n,
\qquad
\gamma_n\in\mathbb{R}.
}
```

The Hilbert–Pólya idea is to explain the reality of `gamma_n` by realizing them as eigenvalues of a self-adjoint operator,

```math
H=H^\ast.
```

The target spectrum is

```math
\operatorname{spec}(H)
=
\{\gamma_1,\gamma_2,\ldots\}.
```

A finite Hermitian matrix fitted to known zeros is not enough. The desired result requires:

- arithmetic construction independent of the target zeros;
- a mathematically defined infinite operator;
- self-adjointness of the limiting operator;
- a rigorous spectral determinant;
- exact equality with the normalized Xi function.

---

# 2. The Xi function

Define

```math
\xi(s)
=
\frac12s(s-1)
\pi^{-s/2}
\Gamma\!\left(\frac{s}{2}\right)
\zeta(s).
```

On the critical line,

```math
\Xi(t)
=
\xi\!\left(\frac12+it\right).
```

For real `t`,

```math
\Xi(-t)=\Xi(t),
```

and the normalized target is

```math
\widehat{\Xi}(t)
=
\frac{\Xi(t)}{\Xi(0)}.
```

The desired identity is

```math
\boxed{
D_H(t)=\widehat{\Xi}(t).
}
```

A finite-grid numerical match is evidence, but not the theorem. The final result must hold as an entire-function identity.

---

# 3. Why the Fredholm determinant is central

For positive spectral values

```math
E_1,E_2,E_3,\ldots,
```

the symmetric spectral determinant formally has the form

```math
D_H(t)
=
\prod_{n=1}^{\infty}
\left(
1-\frac{t^2}{E_n^2}
\right),
```

after the required regularization.

For a finite prefix,

```math
D_N(t)
=
\prod_{n=1}^{N}
\left(
1-\frac{t^2}{E_n^2}
\right).
```

Its zeros occur at

```math
t=\pm E_n.
```

Therefore

```math
D_H(t)
=
\frac{\Xi(t)}{\Xi(0)}
```

would identify the complete Xi zero set with a self-adjoint spectrum. This is much stronger than matching several individual zeros numerically.

---

# 4. Arithmetic source: the von Mangoldt function

The core arithmetic input is the von Mangoldt function,

```math
\Lambda(n)
=
\begin{cases}
\log p, & n=p^k,\\
0, & \text{otherwise}.
\end{cases}
```

It enters naturally through

```math
-\frac{\zeta'(s)}{\zeta(s)}
=
\sum_{n=1}^{\infty}
\frac{\Lambda(n)}{n^s}.
```

The construction therefore retains prime powers, not only bare primes.

For arithmetic cutoff `L`, the frozen density uses weights of the form

```math
w_n(L)
=
\frac{\Lambda(n)}
{\log n\,\sqrt n}
e^{-n/L}.
```

The exponential damping provides a regulated arithmetic problem whose cutoff can later be studied.

---

# 5. Prime phase and Archimedean phase

A representative arithmetic phase is

```math
F_L(E)
=
\frac{\theta(E)}{\pi}
-
\frac{1}{\pi}
\sum_n
w_n(L)\sin(E\log n)
+
C,
```

where `theta(E)` is the Archimedean/Riemann–Siegel phase and `C` fixes the counting convention.

Its derivative is

```math
F_L'(E)
=
\frac{\theta'(E)}{\pi}
-
\frac{1}{\pi}
\sum_n
w_n(L)\log n\cos(E\log n).
```

Substituting the Mangoldt weights gives

```math
F_L'(E)
=
\frac{\theta'(E)}{\pi}
-
\frac{1}{\pi}
\sum_n
\frac{\Lambda(n)}{\sqrt n}
e^{-n/L}
\cos(E\log n).
```

This arithmetic density is constructed without supplying the Riemann zeros.

---

# 6. Quantization law

The arithmetic phase defines consecutive spectral cells. Spectral levels are obtained from a fixed crossing law, schematically

```math
F_L(E_k)=k+\alpha,
```

or equivalently

```math
\Phi_L(E_k)=\pi(k+\beta),
```

depending on the audit convention.

The resulting ordered levels are

```math
E_1^{(L)}
<
E_2^{(L)}
<
\cdots
<
E_N^{(L)}.
```

The important point is that the root solver receives the arithmetic phase law, not a target list of known zeros.

---

# 7. Self-adjoint Jacobi realization

The finite spectral data are converted into a Jacobi operator,

```math
J_N
=
\begin{pmatrix}
a_0 & b_0 & 0 & \cdots \\
b_0 & a_1 & b_1 & \ddots \\
0 & b_1 & a_2 & \ddots \\
\vdots & \ddots & \ddots & \ddots
\end{pmatrix},
```

with

```math
a_n\in\mathbb{R},
\qquad
b_n>0.
```

The baseline operator audit checks:

- quantization residuals;
- self-adjointness defect;
- positivity of Jacobi off-diagonals;
- eigenvalue preservation;
- cutoff stability;
- Dirac/chiral symmetry diagnostics;
- disruption controls.

Because `Xi` is even, the associated Dirac-style symmetric spectrum naturally contains

```math
\pm E_n.
```

The finite even determinant is therefore

```math
D_N(t)
=
\prod_{n=1}^{N}
\left(
1-\frac{t^2}{E_n^2}
\right).
```

---

# 8. ROT interpretation

The ROT component motivates a recursively accumulated positive kernel of the schematic form

```math
K_{m+1}
=
\rho K_m
+
(1-\rho)\,
\kappa_m
|r_m\rangle\langle r_m|.
```

A positive limiting kernel suggests an inverse-square-root Hamiltonian,

```math
H_{\mathrm{ROT}}
=
sK_{\mathrm{ROT}}^{-1/2}.
```

Hence

```math
H_{\mathrm{ROT}}^{-2}
=
\frac1{s^2}K_{\mathrm{ROT}},
```

and the determinant target becomes

```math
\det\!\left(
I-\frac{t^2}{s^2}K_{\mathrm{ROT}}
\right)
=
\frac{\Xi(t)}{\Xi(0)}.
```

Earlier ROT-native realizations were tested and did not automatically satisfy all spectral and convergence requirements. The current project therefore treats ROT as a candidate organizing principle, while requiring independent arithmetic, operator, Fredholm, and falsification tests.

---

# 9. Finite cutoff versus infinite operator

At fixed arithmetic cutoff `L`, the pipeline defines a finite or fixed-`L` spectral operator. The actual Hilbert–Pólya target requires

```math
H_L
\longrightarrow
H_\infty
\qquad
(L\to\infty).
```

There are two distinct limiting problems:

```math
N\to\infty
```

for Jacobi depth, and

```math
L\to\infty
```

for arithmetic cutoff removal.

A complete theorem would require an operator convergence statement such as

```math
(H_L-zI)^{-1}
\longrightarrow
(H_\infty-zI)^{-1}
```

for

```math
z\notin\mathbb{R}.
```

The global cutoff-independent operator theorem is still open.

---

# 10. Why a finite determinant is insufficient

A finite prefix contains only

```math
E_1,\ldots,E_N,
```

whereas `Xi` encodes an infinite spectral set.

Thus the raw determinant

```math
D_N(t)
```

necessarily misses the contribution of spectral states above the last resolved level.

The later ROT–RH pipeline therefore interprets a large part of the finite determinant discrepancy as an unresolved spectral-tail problem rather than immediately retuning the low-energy operator.

This leads to the arithmetic tail completion.

---

# 11. Inverse spectral moments and tail completion

Define

```math
S_m
=
\sum_{n=1}^{\infty}E_n^{-2m},
```

and its finite-prefix version

```math
S_m^{(N)}
=
\sum_{n=1}^{N}E_n^{-2m}.
```

The unresolved contribution is

```math
T_m
=
S_m-S_m^{(N)}.
```

Since

```math
\log D_H(t)
=
-\sum_{m=1}^{\infty}
\frac{S_m}{m}t^{2m},
```

the missing tail contributes

```math
\log D_{\mathrm{tail}}(t)
=
-\sum_{m=1}^{\infty}
\frac{T_m}{m}t^{2m}.
```

At finite moment order `M`, the implemented completed determinant is

```math
\boxed{
D_{\mathrm{completed}}(t)
=
D_N(t)
\exp\!\left[
-\sum_{m=1}^{M}
\frac{T_m}{m}t^{2m}
\right].
}
```

---

# 12. Arithmetic decomposition of the tail

For `p = 2m`, the tail is modeled schematically by

```math
T_m
=
\int_{E_c}^{\infty}
E^{-p}F_L'(E)\,dE
-
\frac12E_c^{-p}.
```

v100.1 decomposes this as

```math
\boxed{
T_m
=
T_m^{\mathrm{Weyl}}
+
T_m^{\mathrm{Arch}}
+
T_m^{\mathrm{prime}}
+
T_m^{\mathrm{mid}}.
}
```

## Weyl term

Using the leading spectral density,

```math
\rho_{\mathrm{Weyl}}(E)
\sim
\frac{1}{2\pi}
\log\!\left(\frac{E}{2\pi}\right),
```

the analytic tail moment is

```math
T_m^{\mathrm{Weyl}}
=
\frac{E_c^{1-2m}}
{2\pi(2m-1)}
\left[
\log\!\left(\frac{E_c}{2\pi}\right)
+
\frac{1}{2m-1}
\right].
```

## Exact Archimedean correction

```math
T_m^{\mathrm{Arch}}
=
\int_{E_c}^{\infty}
\frac{
\theta'(E)/\pi
-
\frac{1}{2\pi}\log(E/2\pi)
}
{E^{2m}}
\,dE.
```

## Prime oscillatory correction

```math
T_m^{\mathrm{prime}}
=
-\frac{1}{\pi}
\sum_n
w_n(L)\log n
\int_{E_c}^{\infty}
\frac{\cos(E\log n)}{E^{2m}}
\,dE.
```

The stable v100.1 implementation evaluates the large-`aE` oscillatory tail using an integration-by-parts asymptotic expansion rather than the earlier unstable upward recurrence.

## Midpoint correction

```math
T_m^{\mathrm{mid}}
=
-\frac12E_c^{-2m}.
```

This half-level term corrects the discrete-to-continuum endpoint.

---

# 13. Frozen ablations

v100.1 constructs four increasingly complete pre-Xi tail laws:

```math
T_m^{(A)}
=
T_m^{\mathrm{Weyl}},
```

```math
T_m^{(B)}
=
T_m^{\mathrm{Weyl}}
+
T_m^{\mathrm{mid}},
```

```math
T_m^{(C)}
=
T_m^{\mathrm{Weyl}}
+
T_m^{\mathrm{Arch}}
+
T_m^{\mathrm{mid}},
```

and

```math
\boxed{
T_m^{(D)}
=
T_m^{\mathrm{Weyl}}
+
T_m^{\mathrm{Arch}}
+
T_m^{\mathrm{prime}}
+
T_m^{\mathrm{mid}}.
}
```

All of these are constructed before Xi is imported.

---

# 14. The Xi firewall

The pipeline is intentionally ordered as

```text
Arithmetic construction
        |
        v
Operator roots
        |
        v
Jacobi / spectral object
        |
        v
Arithmetic tail
        |
        v
Serialize frozen object
        |
        v
SHA-256 hash
        |
   XI FIREWALL
        |
        v
Evaluate Xi
        |
        v
Independent scoring
```

Before the freeze:

- no `Xi(t)` values are used;
- no known zero ordinates are used for parameter selection;
- no tail coefficient is selected from Xi performance.

After the freeze, Xi may be used for determinant, moment, root, or log-derivative scoring.

The manifest records the frozen hash.

---

# 15. Post-freeze Xi moments

The normalized Xi function satisfies locally

```math
\log\!\left(
\frac{\Xi(t)}{\Xi(0)}
\right)
=
-\sum_{m=1}^{\infty}
\frac{S_m^\Xi}{m}t^{2m}.
```

Hence

```math
S_m^\Xi
=
-\frac{m}{(2m)!}
\left.
\frac{d^{2m}}{dt^{2m}}
\log\!\left(
\frac{\Xi(t)}{\Xi(0)}
\right)
\right|_{t=0}.
```

These moments are evaluated only after the arithmetic object is frozen.

---

# 16. Performance

The current repository records the strongest validated values from the v97-v100.1 pipeline.

## Raw determinant

A representative v98 post-freeze run reported

```math
\operatorname{RMS}_{\mathrm{raw}}
\approx
1.564988\times10^{-2}.
```

## Tail-completed determinant

The corresponding reported tail-completed symmetric RMS was

```math
\operatorname{RMS}_{\mathrm{completed}}
\approx
5.793763\times10^{-5}.
```

This is an improvement factor of approximately

```math
\frac{1.564988\times10^{-2}}
{5.793763\times10^{-5}}
\approx
2.70\times10^2.
```

So the determinant discrepancy improved by roughly **270x**.

The important interpretation is that a large fraction of the raw error behaved like missing high-energy spectral mass.

---

# 17. v99 first-moment decomposition

At `L = 16000`, depth 48, v99 reported approximately

```math
E_c
\approx
139.732482236.
```

The finite-prefix first inverse moment was

```math
S_1^{\mathrm{prefix}}
\approx
1.8452179660\times10^{-2}.
```

The leading Weyl tail was

```math
T_1^{\mathrm{Weyl}}
\approx
4.6719998137\times10^{-3}.
```

The Archimedean correction was approximately

```math
T_1^{\mathrm{Arch}}
\approx
-8.1021\times10^{-10},
```

and the prime-density correction was approximately

```math
T_1^{\mathrm{prime}}
\approx
6.39819\times10^{-6}.
```

---

# 18. v100.1 full completion

For depth 32,

```math
E_c
\approx
105.446614,
```

with

```math
S_1^{(32)}
\approx
2.310493074494\times10^{-2}.
```

For depth 48,

```math
E_c
\approx
139.732482,
```

with

```math
S_1^{(48)}
\approx
2.310496887752\times10^{-2}.
```

The depth 32 to 48 full-completion moment convergence was reported as

```math
\frac{
\|S^{(48)}-S^{(32)}\|_2
}{
\|S^{(32)}\|_2
}
\approx
1.650\times10^{-6}.
```

The corresponding completed-determinant convergence was

```math
\frac{
\|D^{(48)}-D^{(32)}\|_2
}{
\|D^{(32)}\|_2
}
\approx
6.442\times10^{-7}.
```

These are **internal convergence diagnostics**. They are not proof-level global errors against Xi.

---

# 19. Performance summary

| Quantity | Reported value | Meaning |
|---|---:|---|
| Raw finite determinant symmetric RMS | `1.564988e-02` | finite-prefix discrepancy |
| Tail-completed representative RMS | `5.793763e-05` | post-completion discrepancy |
| Improvement | `~270x` | missing tail explains substantial raw error |
| Depth 32 -> 48 moment relative L2 | `1.650e-06` | internal moment stability |
| Depth 32 -> 48 determinant relative L2 | `6.442e-07` | internal determinant stability |
| Completed `S1`, depth 32 | `2.310493074494e-02` | full arithmetic completion |
| Completed `S1`, depth 48 | `2.310496887752e-02` | full arithmetic completion |
| Xi used before freeze | `NO` | construction firewall |
| Known zeros used before freeze | `NO` | construction firewall |
| Exact Fredholm-Xi identity proved | `NO` | open theorem |
| RH proved | `NO` | open theorem |

---

# 20. Repository structure

```text
rot-rh-hilbert-polya/
├── .github/
│   └── workflows/
│       └── ci.yml
│
├── docs/
│   ├── THEORY.md
│   ├── OPERATOR_CONSTRUCTION.md
│   ├── FREDHOLM_IDENTITY.md
│   ├── PERFORMANCE.md
│   ├── LIMITATIONS.md
│   ├── ROADMAP.md
│   └── REPRODUCIBILITY.md
│
├── manuscript/
│   ├── ARTICLE.md
│   └── OUTLINE.md
│
├── research/
│   ├── legacy_pipeline/
│   │   ├── rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit.py
│   │   ├── rot_rh_full_rot_hp_pipeline_v97.py
│   │   ├── rot_rh_weyl_tail_completion_v98.py
│   │   ├── rot_rh_tail_midpoint_s1_v99.py
│   │   ├── rot_rh_full_arithmetic_tail_v100_1.py
│   │   └── rot_rh_pnt_finite_part_full_pipeline_v102_1.py
│   │
│   └── red_team/
│       └── rot_rh_adaptive_mesh_surrogates_v129.py
│
├── results/
│   ├── reported_metrics.json
│   └── provenance.json
│
├── scripts/
│   ├── run_smoke.ps1
│   ├── run_smoke.sh
│   ├── run_v97_pipeline.ps1
│   └── run_v100_1.ps1
│
├── src/
│   └── rot_rh_hp/
│       ├── __init__.py
│       ├── arithmetic.py
│       ├── determinant.py
│       ├── metrics.py
│       ├── tail.py
│       └── verification.py
│
├── tests/
│   ├── test_arithmetic.py
│   ├── test_determinant.py
│   ├── test_firewall.py
│   ├── test_manifest.py
│   └── test_tail.py
│
├── CITATION.cff
├── CHANGELOG.md
├── CONTRIBUTING.md
├── LICENSE
├── MANIFEST.sha256
├── pyproject.toml
├── requirements.txt
└── README.md
```

---

# 21. How the repository is organized

The repository deliberately separates four layers.

## Research provenance

`research/legacy_pipeline/` contains the recovered historical scripts. These are retained so the reported numerical results can be traced back to the exact implementations that generated them.

## Reusable mathematical package

`src/rot_rh_hp/` contains smaller reusable implementations of the arithmetic, determinant, tail, metrics, and verification components.

## Tests

`tests/` contains lightweight tests for:

- arithmetic functions;
- determinant identities;
- tail formulas;
- firewall ordering;
- manifest and SHA integrity.

## Documentation

`docs/` separates the theory, operator construction, Fredholm target, limitations, performance, and theorem roadmap.

This separation is intentional: research provenance should not be silently replaced by a cleaner refactor.

---

# 22. Main pipeline versions

## Baseline

```text
rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit.py
```

Builds the Mangoldt prime phase and checks structural Jacobi/Dirac properties.

## v97

```text
rot_rh_full_rot_hp_pipeline_v97.py
```

Main frozen operator pipeline. It constructs and validates pre-Xi gates, freezes the operator, hashes it, and only afterward enables Xi scoring.

## v98

```text
rot_rh_weyl_tail_completion_v98.py
```

Introduces unresolved spectral-tail completion and demonstrates the large determinant improvement.

## v99

```text
rot_rh_tail_midpoint_s1_v99.py
```

Isolates the first inverse spectral moment and verifies the midpoint correction mechanism.

## v100.1

```text
rot_rh_full_arithmetic_tail_v100_1.py
```

Builds the full Weyl + Archimedean + prime + midpoint tail for multiple moments and performs post-freeze Xi comparison.

## v102.1

```text
rot_rh_pnt_finite_part_full_pipeline_v102_1.py
```

Explores a refined PNT finite-part treatment and continuum control.

## v129

```text
rot_rh_adaptive_mesh_surrogates_v129.py
```

Red-team audit using adaptive meshes, permutations, and surrogate arithmetic to detect aliasing or accidental low-frequency structure.

---

# 23. Numerical artifact policy

Earlier exploratory stages uncovered effects that were later traced to interpolation, quadrature, aliasing, or cancellation issues.

The current repository therefore emphasizes:

- precision doubling;
- independent numerical methods;
- grid refinement;
- cutoff refinement;
- depth refinement;
- permutation controls;
- phase scrambling;
- surrogate arithmetic;
- asymptotic checks;
- explicit monitoring of subtractive cancellation.

The goal is not simply to lower RMSE. The goal is to find a mechanism that survives attempts to falsify it.

---

# 24. Installation

Python 3.10+ is recommended.

## Windows PowerShell

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Linux/macOS

```bash
python3.11 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
pip install -r requirements.txt
```

Run the fast tests with

```bash
pytest -q
```

---

# 25. Representative v100.1 command

```powershell
python .\research\legacy_pipeline\rot_rh_full_arithmetic_tail_v100_1.py `
  --v97-prefix .\rot_rh_full_rot_hp_pipeline_v97 `
  --base-module rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit `
  --L 16000 `
  --depths 32,48 `
  --tail-factor 18 `
  --tail-moment-order 10 `
  --post-moment-order 4 `
  --identity-max 50 `
  --identity-step 0.5 `
  --xi-dps 60 `
  --prefix rot_rh_full_arithmetic_tail_v100_1
```

The script verifies the frozen v97 SHA before proceeding.

---

# 26. Exact theorem ladder

The remaining path can be stated as explicit theorem obligations.

## A. Arithmetic phase limit

Prove

```math
F_L(E)\longrightarrow F(E)
```

with adequate uniform control.

## B. Root convergence

Prove

```math
E_n^{(L)}\longrightarrow E_n.
```

## C. Jacobi coefficient convergence

Prove

```math
a_n^{(L)}\to a_n,
\qquad
b_n^{(L)}\to b_n>0.
```

## D. Infinite self-adjoint operator

Construct a dense domain and prove

```math
H_\infty=H_\infty^\ast.
```

## E. Fredholm class

Establish the trace/Schatten properties required to define the regularized determinant of

```math
H_\infty^{-2}.
```

## F. Rigorous arithmetic-tail error bound

Prove a bound of the form

```math
\left|
T_m^{\mathrm{discrete}}
-
T_m^{\mathrm{model}}
\right|
\le
\varepsilon_m(L,E_c),
```

with

```math
\varepsilon_m(L,E_c)\to0.
```

## G. Compact determinant convergence

Prove

```math
D_L(t)\longrightarrow D_\infty(t)
```

uniformly on compact subsets of the complex plane.

## H. Exact Xi identity

Finally prove

```math
\boxed{
D_\infty(t)
=
\frac{\Xi(t)}{\Xi(0)}.
}
```

Only after these steps would the Hilbert–Pólya route be complete.

---

# 27. Operator/Fredholm formulation through the inverse kernel

Define

```math
K_\infty=H_\infty^{-2}.
```

Then

```math
\operatorname{spec}(K_\infty)
=
\left\{
E_1^{-2},
E_2^{-2},
\ldots
\right\},
```

and

```math
D_\infty(t)
=
\det\!\left(
I-t^2K_\infty
\right).
```

The inverse moments become traces:

```math
S_m
=
\operatorname{Tr}(K_\infty^m).
```

Therefore

```math
\log\det(I-t^2K_\infty)
=
-\sum_{m=1}^{\infty}
\frac{t^{2m}}{m}
\operatorname{Tr}(K_\infty^m).
```

This is the clean mathematical bridge between the operator and the Xi moment hierarchy.

---

# 28. What the repository currently supports

The strongest defensible claims are:

1. A zero-free Gamma + Mangoldt arithmetic phase can generate a structurally valid self-adjoint finite Jacobi spectrum.
2. The raw finite determinant contains a large missing-tail error.
3. Much of that error is explained by independently derived high-energy spectral density.
4. The unresolved tail decomposes naturally into Weyl, Archimedean, prime, and midpoint pieces.
5. The completed determinant shows strong depth stability in the tested regime.
6. The code enforces a pre-freeze separation between construction and Xi scoring.
7. The remaining problem has been reduced to explicit operator-limit, tail-bound, and entire-function identity theorems.

---

# 29. What the repository does not claim

The repository does not claim that:

- a low RMSE proves RH;
- the finite Jacobi matrix is already the final infinite operator;
- the arithmetic cutoff limit has been proved;
- the tail approximation is already rigorously uniform in the required joint limit;
- the completed determinant has been proved to converge on all complex compact sets;
- the exact Xi identity has been established.

The current theorem status is therefore

```math
\boxed{
\text{Riemann Hypothesis proof status: OPEN}.
}
```

---

# 30. Final perspective

The strongest current formulation of the program is

```math
\boxed{
\Lambda(n)
\longrightarrow
F_L(E)
\longrightarrow
\{E_n^{(L)}\}
\longrightarrow
H_L
\longrightarrow
D_L(t)
\longrightarrow
D_L^{\mathrm{completed}}(t)
\longrightarrow
\frac{\Xi(t)}{\Xi(0)}.
}
```

Every arrow through the completed determinant can now be numerically tested.

The final arrow is still a numerical comparison rather than a theorem.

The central open question is

```math
\boxed{
\text{Does the zero-free arithmetic sequence converge to a self-adjoint }
H_\infty
\text{ whose exact Fredholm determinant is }
\Xi(t)/\Xi(0)?
}
```

At present,

```math
\boxed{
\text{strong numerical evidence}
\neq
\text{proof}.
}
```

The purpose of this repository is to turn that remaining gap into a precise, reproducible theorem-development program.

---

# 31. Current implementation status

```text
Finite arithmetic operator construction       : IMPLEMENTED
Finite self-adjoint Jacobi realization         : IMPLEMENTED
Zero-free pre-freeze construction              : IMPLEMENTED
SHA-256 freeze firewall                        : IMPLEMENTED
Weyl tail completion                           : IMPLEMENTED
Archimedean tail correction                    : IMPLEMENTED
Prime oscillatory tail correction              : IMPLEMENTED
Midpoint correction                            : IMPLEMENTED
Multi-moment Fredholm completion               : IMPLEMENTED
Post-freeze Xi scoring                         : IMPLEMENTED
Depth convergence audits                       : IMPLEMENTED
Adaptive/surrogate falsification tests         : IMPLEMENTED

Global arithmetic cutoff limit                 : OPEN
Infinite self-adjoint operator theorem         : OPEN
Rigorous arithmetic-tail remainder bound       : OPEN
Locally uniform complex Fredholm convergence   : OPEN
Exact Fredholm-Xi identity                     : OPEN
Riemann Hypothesis proof                       : OPEN
```
