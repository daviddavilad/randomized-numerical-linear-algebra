# Predicting q=1 from the q=0 curve

*Derived 2026-08-22, before adding q to the benchmark.*

## Claim

Power iteration raises $\sigma_i \mapsto \sigma_i^{2q+1}$, so the gap at the
truncation point becomes $g^{2q+1}$ where $g = \sigma_{k+1}/\sigma_k$. If the
$q=0$ error curve were a function of the gap alone, then the $q=1$ ratio for a
spectrum with gap $g$ should equal the $q=0$ ratio of a spectrum with gap $g^3$.

## Lookup

$q=0$ measurements, $m=300$, $n=200$, $k=20$, $p=10$, 10 seeds:

| spectrum | gap | mean ratio |
|---|---|---|
| exp 1.0 | 0.3679 | 1.000000105 |
| exp 0.5 | 0.6065 | 1.000530152 |
| exp 0.3 | 0.7408 | 1.013998225 |
| exp 0.2 | 0.8187 | 1.062259422 |
| poly 3.0 | 0.8638 | 1.227023766 |
| exp 0.1 | 0.9048 | 1.198479169 |
| poly 2.0 | 0.9070 | 1.318590745 |
| poly 1.5 | 0.9294 | 1.337235016 |
| poly 1.0 | 0.9524 | 1.298280110 |
| poly 0.5 | 0.9759 | 1.163490324 |
| flat | 1.0000 | 1.000000000 |

Predicted $q=1$ by looking up $g^3$:

| spectrum | $g$ | $g^3$ | nearest row | predicted $q=1$ |
|---|---|---|---|---|
| exp 0.1 | 0.9048 | 0.7407 | exp 0.3 | $\approx 1.014$ |
| poly 2.0 | 0.9070 | 0.7462 | exp 0.3 | $\approx 1.014$ |
| poly 1.5 | 0.9294 | 0.8028 | between exp 0.3 and exp 0.2 | $\approx 1.04$ |
| poly 1.0 | 0.9524 | 0.8637 | poly 3.0 | $\approx 1.23$ (?) |
| poly 0.5 | 0.9759 | 0.9294 | poly 1.5 | $\approx 1.34$ (?) |

## Where the lookup breaks

The last two rows are inconsistent with a theorem. Monotonicity in $q$ was
verified in `power-iteration-stability.md`; poly 0.5 measures $1.163$ at $q=0$,
so its $q=1$ value cannot be $1.34$.

The lookup is thus invalid, and the reason is the result from
`derivations/notes/` on the gap hypothesis: the $q=0$ ratio is **not** a
function of the gap. Two spectra with matching gaps gave $1.198$ and $1.319$
(exp 0.1 vs poly 2.0, gaps $0.9048$ and $0.9070$). Power iteration changes
$g \mapsto g^3$ but it also changes the whole tail - the tail dimension
$\mathrm{sr}_{\text{tail}} = \left(\sum_{j>k}\sigma_j^2\right)/\sigma_{k+1}^2$
is computed from $\sigma^{2q+1}$, and a flatter spectrum's tail contracts far
less than a steeply decaying one's.

As a result, the lookup should hold where the gap dominates and fail where
the tail does.

## Falsifiable version

- Rows with small $g$ (exp 0.3 and below): lookup accurate, all $\to 1.00$.
- Rows in the middle (poly 2.0, poly 1.5): lookup roughly accurate.
- Rows with $g \to 1$ (poly 0.5, flat): lookup fails; measured values must
  respect monotonicity, so they will be **below** the lookup prediction.

If poly 0.5 comes in near $1.16$ or lower rather than $1.34$, the tail
explanation is correct and the gap-only claim is falsified.

## Outcome

Measured $q=1$ against the predictions above:

| spectrum | $g$ | predicted | measured |
|---|---|---|---|
| exp 0.1 | 0.9048 | $\approx 1.014$ | 1.000519 |
| poly 2.0 | 0.9070 | $\approx 1.014$ | 1.001151 |
| poly 1.5 | 0.9294 | $\approx 1.04$ | 1.003021 |
| poly 1.0 | 0.9524 | $\approx 1.23$ (?) | 1.007238 |
| poly 0.5 | 0.9759 | $\approx 1.34$ (?) | 1.015473 |

The two rows flagged as inconsistent with monotonicity are confirmed invalid:
poly 0.5 measures $1.0155$, not $1.34$. But the lookup also fails on the rows
it was expected to get right - poly 2.0 came in $12\times$ better than
predicted.

Both failures have an explanation. Amplification contracts the tail as well as
widening the gap: at $q=1$, poly $\alpha=2.0$ becomes $\sigma_i^3 = (i+1)^{-6}$
with $\mathrm{sr}_{\text{tail}} \approx 1.08$, while exp $\alpha=0.3$ - the
lookup match on gap - has $\mathrm{sr}_{\text{tail}} \approx 2.2$. Matching on
gap alone compares against a strictly harder problem, so the lookup is
*underestimating* the benefit wherever the tail collapses, and overestimating it
where the spectrum is too flat for cubing to move anything.

## Power iteration restores gap-monotonicity

Excess over optimal, $\varepsilon_q = \text{ratio} - 1$, sorted by gap:

| gap | $\varepsilon_0$ | $\varepsilon_1$ |
|---|---|---|
| 0.8638 (poly 3.0) | 2.27e-1 | 1.49e-4 |
| 0.9048 (exp 0.1) | **1.99e-1** | 5.20e-4 |
| 0.9070 (poly 2.0) | 3.19e-1 | 1.15e-3 |
| 0.9294 (poly 1.5) | 3.37e-1 | 3.02e-3 |
| 0.9524 (poly 1.0) | **2.98e-1** | 7.24e-3 |
| 0.9759 (poly 0.5) | **1.64e-1** | 1.55e-2 |

At $q=0$ the ordering is broken twice: poly 3.0 sits above exp 0.1 despite a
smaller gap, and the curve peaks at $0.9294$ and then declines. At $q=1$ it is
monotone increasing across the whole range, with the two families interleaving
in gap order and no inversions. Same at $q=2$.

The non-monotonic structure that falsified the local-gap hypothesis is therefore
simply a $q=0$ case. One power iteration makes the gap the governing variable it
was not before.

## Measurement limit

exp $\alpha=1.0$ at $q=2$ reads $0.999999999$, which is below the Eckart–Young
optimum, *impossible* by construction. The optimum is computed
analytically from the known spectrum, but the achieved error is computed
numerically, and at the $10^{-9}$ level they cross. Values below $\sim 10^{-8}$ in
these columns are measurement noise.