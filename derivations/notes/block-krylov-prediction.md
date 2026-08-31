# Block Krylov - gap-independence prediction

*Written 2026-08-29, before implementing block Krylov.*

## The two methods

Power iteration keeps only the final iterate of $(AA^\top)^q A\Omega$. Block
Krylov keeps every one:

$$K = \big[\,A\Omega,\; (AA^\top)A\Omega,\; \ldots,\; (AA^\top)^q A\Omega\,\big]$$

Same passes over $A$ at the same $q$. Different subspace: $(q{+}1)\ell$ columns
instead of $\ell$, so block Krylov buys accuracy with memory rather than passes.

## Claim

Musco & Musco (2015) give block Krylov an $O(1/\sqrt{\varepsilon})$ iteration
count against power iteration's $O(1/\varepsilon)$, and, the part that matters
here, the bound does not depend on the spectral gap.

## Prediction

If gap-independence truly exists and is not just an artifact of the analysis, it should
show up as a **flatter error curve**, not a uniform improvement.

Power iteration at $q=1$, excess over optimal:

| gap | $\varepsilon_1$ |
|---|---|
| 0.3679 | 1.0e-9 |
| 0.8638 | 1.49e-4 |
| 0.9070 | 1.15e-3 |
| 0.9759 | 1.55e-2 |

Seven orders of magnitude across the sweep. **Block Krylov at $q=1$ should
compress that range.** A uniform downward shift (every row improving by
roughly the same factor) would mean the gap still governs and the
gap-independence is not visible at these sizes.

Secondary: the flat spectrum should stay at exactly 1.000000000, since the
exactness argument does not depend on which subspace is used, only that it lies
in $\mathrm{range}(A)$.

## Outcome

*To be filled after measurement.*