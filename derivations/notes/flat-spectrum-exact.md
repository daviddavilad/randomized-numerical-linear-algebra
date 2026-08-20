# Flat spectrum - a special case for randomized SVD

*Derived 2026-08-20, before running the flat benchmark case. For measurements, 
see `benchmarks/bench_spectrum.cpp`.*

A rare case where the rSVD error can be computed in closed form, with no
probabilistic bound. Used as a falsification test for the local-gap hypothesis.

## Setup

Let $A \in \mathbb{R}^{m \times n}$ with $r = \min(m,n)$ and a **flat**
spectrum, $\sigma_1 = \cdots = \sigma_r = 1$:

$$A = U V^\top, \qquad U^\top U = V^\top V = I_r .$$

Then $\|A\|_F^2 = r$, and the Eckart–Young optimum at rank $k$ is

$$\|A - A_k\|_F = \Big(\sum_{i>k} \sigma_i^2\Big)^{1/2} = \sqrt{r-k}.$$

The algorithm, with sketch width $\ell = k+p$:

$$Y = A\Omega, \quad Q = \operatorname{orth}(Y), \quad B = Q^\top A,
\quad \widehat{A} = Q B_k$$

where $B_k$ is the rank-$k$ truncated SVD of $B$. Assume $\ell \le r$ and
$\operatorname{rank}(Y) = \ell$ (true with probability 1).

## Step 1 - split the error

$$A - QB_k = \underbrace{(I - QQ^\top)A}_{\perp\,\mathrm{range}(Q)} + \underbrace{Q(B - B_k)}_{\in\,\mathrm{range}(Q)}$$

The two terms have columns in orthogonal subspaces, so

$$\|A - QB_k\|_F^2 = \|(I-QQ^\top)A\|_F^2 + \|B - B_k\|_F^2 ,$$

using $\|QX\|_F = \|X\|_F$ for $Q$ with orthonormal columns.

## Step 2 - range of $Q$

$Y = A\Omega$, so every column of $Y$ is a combination of columns of $A$:

$$\mathrm{range}(Q) = \mathrm{range}(Y) \subseteq \mathrm{range}(A)
= \mathrm{range}(U) \quad\Longrightarrow\quad UU^\top Q = Q .$$

## Step 3 - $B$ has $\ell$ singular values, all equal to 1

$$BB^\top = Q^\top A A^\top Q = Q^\top U \Sigma^2 U^\top Q = Q^\top U U^\top Q = Q^\top Q = I_\ell ,$$

using $\Sigma = I_r$ and Step 2. So $\sigma_i(B) = 1$ for $i = 1,\dots,\ell$,
and

$$\|B\|_F^2 = \ell, \qquad \|B - B_k\|_F^2 = \ell - k = p .$$

## Step 4 - the projection term

$$\|(I-QQ^\top)A\|_F^2 = \|A\|_F^2 - \|QQ^\top A\|_F^2
 = r - \|Q^\top A\|_F^2 = r - \|B\|_F^2 = r - \ell .$$

## Result

$$\|A - QB_k\|_F^2 = (r-\ell) + (\ell-k) = r-k$$

$$\boxed{\;\frac{\|A - QB_k\|_F}{\|A - A_k\|_F}
= \frac{\sqrt{r-k}}{\sqrt{r-k}} = 1\;}$$

**Randomized SVD is exactly optimal on a flat spectrum**, for every draw of
$\Omega$ - not approximately, not in expectation. The oversampling $p$ cancels:
the $p$ extra sketch directions are exactly the $p$ discarded in truncation.

## Reason

Every $k$-dimensional subspace of $\mathrm{range}(A)$ captures exactly $k$
units of energy, so there is no "correct" top-$k$ subspace to miss. The
randomness has nothing to get wrong.

## Prediction for the benchmark

| hypothesis | quantity | prediction |
|---|---|---|
| local gap | $\sigma_{k+1}/\sigma_k = 1.000$ | worst ratio in the study, $> 1.34$ |
| tail dimension | $\mathrm{sr}_{\text{tail}} = r-k = 180$ | ratio $\approx 1.000$ |

The derivation above gives $1.000$ exactly, and `mean` and `max` agree to
machine precision across all seeds.

## Outcome

Measured at $m=300$, $n=200$, $k=20$, $p=10$, over 10 seeds:

| spectrum | gap | mean | max |
|---|---|---|---|
| flat | 1.0000 | 1.000000000 | 1.000000000 |
| exp $\alpha=1.0$ | 0.3679 | 1.000000105 | 1.000000675 |

The flat row shows **zero variance across all ten seeds**, which is the
observable signature of a draw-independent result - the derivation says the
error does not depend on $\Omega$, and the measurement agrees to nine decimals.

The contrast with exp $\alpha=1.0$ is the point. Both print `1.0000` at four
decimals, but that row has visible seed-to-seed scatter: it is *approximately*
optimal because its tail is tiny ($\|A-A_k\|_F/\|A\|_F \approx 2\times10^{-9}$),
not *exactly* optimal.