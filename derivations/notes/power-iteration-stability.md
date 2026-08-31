# Power iteration - gap amplification and its numerical cost

*Derived 2026-08-21. The crossover table was worked out from the spectrum
before re-running the benchmark.*

## Gap amplification

Power iteration replaces $Y = A\Omega$ with $Y = (AA^\top)^q A\,\Omega$.
Writing $A = U\Sigma V^\top$,

$$(AA^\top)^q A = U\Sigma^{2q+1}V^\top ,$$

so the range is unchanged while every singular value is raised to $2q+1$. The
gap the range finder must resolve becomes

$$\left(\frac{\sigma_{k+1}}{\sigma_k}\right)^{2q+1}.$$

## The exponent widens the dynamic range

Across the $\ell = k+p$ retained directions,

$$\frac{\sigma_1^{2q+1}}{\sigma_\ell^{2q+1}}
= \left(\frac{\sigma_1}{\sigma_\ell}\right)^{2q+1}.$$

Once this exceeds $1/\varepsilon \approx 4.5\times10^{15}$, the smallest
directions of $Y$ lie below the rounding error of the largest: $Y$ is
numerically rank-deficient though mathematically full rank, and
$\operatorname{orth}(Y)$ returns meaningless columns for the lost directions.

For $\sigma_i = (i+1)^{-2}$ with $\ell = 30$, the base ratio is $30^2 = 900$:

| $q$ | $900^{2q+1}$ | margin to $1/\varepsilon$ |
|---|---|---|
| 0 | $9.0\times10^{2}$ | $10^{13}$ |
| 1 | $7.3\times10^{8}$ | $10^{7}$ |
| 2 | $5.9\times10^{14}$ | $8\times$ |
| 3 | $4.8\times10^{20}$ | exceeded by $10^{5}$ |

## Fix

Re-orthonormalize after every application of $A^\top$ and of $A$
(HMT Alg. 4.4):

$$Y \leftarrow \operatorname{orth}(A^\top Y), \qquad
  Y \leftarrow \operatorname{orth}(A\,Y).$$

The subspace is unchanged, only the representation is reconditioned, resetting
the dynamic range to $O(1)$ at each half-step so it never accumulates. Both
calls are required - the collapse can occur within one application.

## Measured

poly $\alpha=2.0$, $m=300$, $n=200$, $k=20$, $p=10$, seed 42:

| $q$ | naive | re-orthonormalized |
|---|---|---|
| 0 | 1.273645689 | 1.273645689 |
| 1 | 1.000624473 | 1.000624473 |
| 2 | 1.000012155 | 1.000004883 |
| 3 | 1.645690087 | 1.000000050 |

Consistent with the table: identical where the margin is large, $2.5\times$
degraded at $q=2$ where the margin is $8$, divergent at $q=3$.

## Cross-check

The stable values follow from the amplified gap alone. For $\alpha=2.0$,
$\sigma_{k+1}/\sigma_k = 0.907$, giving $0.907^3 = 0.746$ at $q=1$ and
$0.907^5 = 0.614$ at $q=2$. The $q=0$ sweep at comparable gaps ($0.741$,
$0.607$) measured $1.014$ and $1.0005$; power iteration measures $1.00062$ and
$1.0000049$.