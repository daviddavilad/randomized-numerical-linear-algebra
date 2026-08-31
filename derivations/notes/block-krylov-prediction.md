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

Mean excess over optimal, 10 seeds, $`m=300`$, $`n=200`$, $`k=20`$, $`p=10`$.
Both methods at $`q=1`$ - three passes over $`A`$ either way.

| spectrum | gap | $`\mathrm{sr}_\mathrm{tail}`$ | power iter. | block Krylov | factor |
|---|---|---|---|---|---|
| poly 3.0 | 0.8638 | 4.72 | 1.49e-4 | 3.70e-6 | 40× |
| exp 0.1 | 0.9048 | 5.52 | 5.20e-4 | 1.90e-6 | 274× |
| poly 2.0 | 0.9070 | 7.51 | 1.15e-3 | 1.34e-4 | 8.6× |
| poly 1.5 | 0.9294 | 10.90 | 3.02e-3 | 6.68e-4 | 4.5× |
| poly 1.0 | 0.9524 | 19.31 | 7.24e-3 | 2.61e-3 | 2.8× |
| poly 0.5 | 0.9759 | 47.89 | 1.55e-2 | 6.32e-3 | 2.4× |

Block Krylov wins on every row. **The prediction is rejected**: across the
polynomial rows the power-iteration excess spans 104×, the block Krylov excess
spans 1708×. The curve got *steeper*, not flatter, and the improvement shrinks
monotonically as the problem gets harder.

The improvement factor is not a function of the gap or of
$`\mathrm{sr}_\mathrm{tail}`$: exp 0.1 ($`\mathrm{sr}_\mathrm{tail} = 5.52`$)
improves 274× while poly 3.0 ($`\mathrm{sr}_\mathrm{tail} = 4.72`$) improves 40×.
Same shape of family anomaly that rejected H1.

## Prediction error

The mistake was in the prediction, not in the theorem. Musco & Musco's
gap-independence is a statement about **iteration count to reach a target
$`\varepsilon`$** - block Krylov needs $`q = O(1/\sqrt{\varepsilon})`$ where
power iteration needs $`O(1/\varepsilon)`$, and neither bound involves the gap.
It says nothing about error at *fixed* $`q`$, which is what was measured here.

Testing the actual claim requires fixing a target error and measuring the $`q`$
each method needs to reach it. Follow-up item for the future.