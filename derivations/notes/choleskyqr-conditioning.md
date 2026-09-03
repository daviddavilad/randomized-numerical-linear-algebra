# CholeskyQR and the conditioning threshold

*Written 2026-09-02, before running the sweep.*

## The mechanism

CholeskyQR computes $`G = A^\top A`$, factors $`G = R^\top R`$, and solves $`Q = AR^{-1}`$. Three BLAS-3 operations against Householder's sequential reflector sweep.

The Gram matrix squares the spectrum: $`\sigma_i(G) = \sigma_i(A)^2`$, so $`\kappa(G) = \kappa(A)^2`$. The Cholesky therefore sees a problem quadratically harder than the one posed.

## Threshold

A Cholesky factorization loses positive definiteness once $`\kappa(G)`$ approaches $`1/\varepsilon \approx 4.5\times10^{15}`$. In terms of the input,

$$\kappa(A) \gtrsim \varepsilon^{-1/2} \approx 6.7\times10^{7}.$$

## Predictions

1. Householder stays at $`\sim10^{-15}`$ for every $`\kappa`$ tested - it is
   backward stable and never forms $`G`$.
2. CholeskyQR matches it while $`\kappa \ll 10^8`$.
3. Orthogonality degrades **before** outright failure, as with power iteration
   at $`q=2`$: expect visible loss around $`\kappa \sim 10^6`$–$`10^7`$, where
   $`\kappa^2`$ is within a few orders of $`1/\varepsilon`$.
4. `dpotrf` throws somewhere near $`\kappa \sim 10^8`$.

Rough expectation, orthogonality error against $`\kappa`$:

| $`\kappa`$ | $`\kappa^2`$ | Householder | CholeskyQR |
|---|---|---|---|
| 1e2 | 1e4 | 1e-15 | 1e-15 |
| 1e4 | 1e8 | 1e-15 | 1e-15 |
| 1e6 | 1e12 | 1e-15 | ~1e-12 |
| 1e8 | 1e16 | 1e-15 | fails or $`O(1)`$ |
| 1e10 | 1e20 | 1e-15 | fails |

## Outcome

$`m = 200`$, $`n = 30`$, exponential spectra with $`\alpha = \ln\kappa/(n-1)`$. Achieved conditioning matched the target exactly in every row.

| $`\kappa`$ | $`\kappa^2`$ | Householder | CholeskyQR |
|---|---|---|---|
| 1e2 | 1e4 | 2.63e-15 | 6.52e-13 |
| 1e4 | 1e8 | 4.13e-15 | 2.54e-09 |
| 1e6 | 1e12 | 3.49e-15 | 2.06e-05 |
| 1e7 | 1e14 | 4.09e-15 | 1.07e-03 |
| 1e8 | 1e16 | 2.68e-15 | 1.52e-01 |
| 1e10 | 1e20 | 3.20e-15 | fails, minor 25 |
| 1e12 | 1e24 | 3.61e-15 | fails, minor 22 |

Prediction 1 holds: Householder is flat at $`\sim3\times10^{-15}`$ across ten orders of magnitude of conditioning.

**Predictions 2 and 3 were wrong.** There is no regime where CholeskyQR matches Householder, and no threshold at which degradation begins. Dividing through by $`\varepsilon\kappa^2`$:

| $`\kappa`$ | $`\varepsilon\kappa^2`$ | measured | ratio |
|---|---|---|---|
| 1e2 | 2.2e-12 | 6.5e-13 | 0.29 |
| 1e4 | 2.2e-8 | 2.5e-9 | 0.11 |
| 1e6 | 2.2e-4 | 2.1e-5 | 0.093 |
| 1e7 | 2.2e-2 | 1.1e-3 | 0.048 |
| 1e8 | 2.2 | 1.5e-1 | 0.068 |

Constant within a factor of six over six decades. The law is

$$\|Q^\top Q - I\|_F \approx c\,\varepsilon\,\kappa(A)^2, \qquad c \approx 0.1,$$

holding from the smallest $`\kappa`$ tested. The error is quadratic in $`\kappa`$ *everywhere* - the $`\varepsilon^{-1/2}`$ figure is not where degradation starts but where $`\varepsilon\kappa^2`$ reaches $`O(1)`$ and the accumulated error has consumed the entire answer.

Prediction 4 is also slightly off: at $`\kappa = 10^8`$ the factorization does not throw - it returns $`1.5\times10^{-1}`$, i.e. garbage without complaint. The throw arrives at $`10^{10}`$. Silent failure precedes loud failure, as with naive power iteration at $`q=2`$.

The leading minor at which `dpotrf` breaks down (25 at $`\kappa=10^{10}`$, 22 at $`10^{12}`$) moves earlier with worse conditioning, giving a rough proxy for how far past the edge the input sits.

## Next: CholeskyQR2

Running the algorithm twice - $`Q_1 = \mathrm{cholqr}(A)`$, then $`Q = \mathrm{cholqr}(Q_1)`$ - is claimed to restore machine-precision orthogonality for $`\kappa(A)`$ up to roughly $`\varepsilon^{-1/2}`$.

The mechanism to verify before measuring: $`Q_1`$ is not orthonormal, but its *conditioning* is far better than $`A`$'s. If $`\kappa(Q_1) = O(1)`$ whenever the first pass completes at all, the second pass is a well-conditioned problem and returns machine precision. The open question is whether $`\kappa(Q_1)`$ tracks.