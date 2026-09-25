
# Sparse Sketches: Accuracy, Sparsity, and Rotation

## Objective

Investigate how the sparsity parameter $\zeta$ affects the accuracy of sparse sign sketches, and whether changing the orientation of the input matrix changes the observed results.

The investigation covers three questions:

1. Does increasing $\zeta$ reduce the global Gram-matrix distortion?
2. Is there a measurable accuracy penalty for sparse sign sketches relative to Gaussian sketches?
3. Does rotating the input matrix change the comparison?

The results concern **approximation accuracy**, not computational performance.

---

## 1. The $\zeta$-Independence Derivation

### Setup

Let

$$
S\in\mathbb R^{m\times n}
$$

be a sparse sign sketch with exactly $\zeta$ nonzero entries per column.

Each nonzero entry has magnitude

$$
\frac{1}{\sqrt{\zeta}}
$$

and an independent random sign.

Consequently, every column has unit Euclidean norm:

$$
\|S(:,j)\|_2^2=1.
$$

We study the Gram matrix

$$
G=S^TS
$$

and its deviation from the identity:

$$
G-I.
$$

### Derivation

For two distinct columns $j\neq k$, let $r$ denote the number of overlapping nonzero coordinates.

The expected overlap is

$$
\mathbb E[r]=\frac{\zeta^2}{m}.
$$

Conditioned on $r$, the inner product is a sum of $r$ random signed terms:

$$
G_{jk}
=
\frac{1}{\zeta}
\sum_{t=1}^{r}\varepsilon_t,
$$

where the $\varepsilon_t$ are independent Rademacher random variables.

Therefore,

$$
\mathbb E[G_{jk}\mid r]=0
$$

and

$$
\mathrm{Var}(G_{jk}\mid r)
=
\frac{r}{\zeta^2}.
$$

Taking expectations over the overlap distribution gives

$$
\begin{aligned}
\mathbb E[G_{jk}^2]
&=
\mathbb E\left[\frac{r}{\zeta^2}\right]\\
&=
\frac{1}{\zeta^2}
\frac{\zeta^2}{m}\\
&=
\boxed{\frac{1}{m}}.
\end{aligned}
$$

**The expected squared off-diagonal inner product is independent of $\zeta$.**

Because the columns are exactly normalized, the diagonal entries satisfy

$$
G_{jj}=1.
$$

Hence,

$$
\begin{aligned}
\mathbb E\left[\|S^TS-I\|_F^2\right]
&=
\sum_{j\neq k}\mathbb E[G_{jk}^2]\\
&=
\boxed{\frac{n(n-1)}{m}}.
\end{aligned}
$$

The resulting root-mean-square scale is

$$
\boxed{
\sqrt{\mathbb E\left[\|S^TS-I\|_F^2\right]}
=
\sqrt{\frac{n(n-1)}{m}}.
}
$$

This is a prediction for the square root of the expected **squared** Frobenius error, not an exact identity for $\mathbb E[\|S^TS-I\|_F]$.

### Prediction

Increasing $\zeta$ should not systematically reduce this global Frobenius distortion measure.

Greater sparsity changes the distribution of individual inner products, but the increase in expected overlap is exactly offset by the smaller magnitude of each overlapping contribution.

### Measurement

Parameters: $m=200$, $n=30$, $\zeta=8$ for the reported Frobenius distortion; the sparsity sweep used the same $m$ and $n$.

| Quantity | Result |
|---|---:|
| Predicted RMS scale | 2.086 |
| Measured Frobenius distortion | 2.039 |
| Tested sparsities | $\zeta\in\{2,\ldots,32\}$ |

The measured distortion was approximately flat across the tested values of $\zeta$ and matched a unit-normalized Gaussian sketch to within 2% in a single comparison, with no error bars on either.

The agreement supports the second-moment prediction.

It does **not** imply that sparse sign and Gaussian sketches have identical distributions, tail behavior, spectral errors, or performance on every downstream task.

---

## 2. H8: Does Sparsity Introduce an Accuracy Penalty?

### Hypothesis

H8 investigated whether sparse sign sketches incur a measurable accuracy penalty relative to the Gaussian reference, and whether that penalty decreases as $\zeta$ increases.

### Parameters

$m=300$, $n=200$, $k=20$, $p=10$, $q=0$, across four spectra (exponential $\alpha\in\{0.3,0.1\}$, polynomial $\alpha\in\{2.0,1.0\}$). Errors are measured as the ratio to the Eckart–Young optimum.

### Initial measurement: 20 seeds

At 20 seeds, 15 of 16 comparisons had a positive sign, suggesting a small penalty.

This was preliminary evidence, not a stable conclusion.

### Expanded measurement: 200 seeds

After increasing the experiment to 200 seeds:

- $\zeta=1$ remained approximately $8$–$12$ standard errors worse.
- Every tested $\zeta\geq2$ was within $1.3$ standard errors of the reference.
- The sign pattern reversed: 12 of 16 comparisons were negative.

The reversal matters. A consistent effect should generally become easier to resolve as the sample size increases, rather than reversing direction when the estimate becomes more precise.

### Conclusion

The experiment detected a clear accuracy penalty at

$$
\boxed{\zeta=1}.
$$

For

$$
\boxed{\zeta\geq2},
$$

no systematic accuracy penalty was resolved at 200 seeds.

This is a **failure to detect a penalty**, not proof of equivalence between sparse sign and Gaussian sketches. A formal equivalence claim would require a prespecified tolerance and an appropriate equivalence analysis.

The 20-seed result illustrates why small apparent effects should not be treated as established findings before checking their stability.

---

## 3. Rotation Experiment

### Motivation

Sparse sign sketches are not generally rotation-invariant.

Therefore, results obtained with one fixed orientation of the input matrix need not describe every possible alignment between the input subspace and the sketch coordinates.

The original fixed-rotation experiment tested only one alignment.

### Experimental change

The experiment was repeated while varying the input rotation.

The comparison retained a paired design so that the effect of changing the rotation could be assessed through differences between the two sketching methods.

### Prediction

If orientation materially affects the relative accuracy of the methods, changing the rotation could alter the estimated paired difference or reveal an effect hidden by the original alignment.

### Measurement

After varying the rotation:

- 15 of 16 standard errors decreased.
- No rotation effect appeared in the means.
- The extra variance an alignment effect would be expected to produce in the sparse columns did not appear either.

### Interpretation

Pairing was already in the design before the rotation varied, and it produced no reduction then: with the matrix fixed, the two methods shared no source of variation for the difference to cancel. Once the rotation varied and was shared between the two methods within a seed, there was a genuine common component, and the standard errors fell.

This is the same point as the seed discussion in Section 4, arriving from the other direction. Pairing helps exactly when a real source of variation is shared.

The experiment did not identify a measurable orientation-dependent difference at the tested sample size.

This is a **null result for the tested design**, not a general proof of rotation invariance. Sparse sign sketches remain non-rotation-invariant as a distribution.

---

## 4. Scope and Limitations

### Accuracy, not computational cost

These experiments measure approximation accuracy.

They do not establish the computational advantage of sparse sign sketches because the implementation returns the sketch as a dense matrix.

Consequently, increasing or decreasing $\zeta$ does not currently provide a valid measurement of the computational benefits of sparse storage and sparse matrix operations.

A separate benchmark is needed to measure:

- Sketch construction time.
- Matrix multiplication time.
- Memory usage.
- Scaling with $m$, $n$, and $\zeta$.

The implementation would need to preserve sparse structure rather than immediately materializing a dense matrix.

### Pairing and common randomness

Pairing reduced variance when both methods shared the same randomized input rotation.

However, when the matrix was fixed, using the same integer seed did not produce the same benefit.

The reason is that the two sketch generators use different random transformations. Passing the same seed integer to both generators does not ensure that they share the same underlying random variables.

Thus,

$$
\boxed{
\text{same seed}\neq\text{shared randomness}.
}
$$

Effective pairing requires an actual shared source of experimental variation, not merely identical seed labels.

### Limits of the $\zeta$-independence result

The derivation establishes

$$
\mathbb E\left[\|S^TS-I\|_F^2\right]
=
\frac{n(n-1)}{m}
$$

under the specified sparse sign construction.

It does not establish $`\zeta`$-independence for every accuracy metric.

In particular, the result does not determine:

- Spectral-norm distortion.
- Tail probabilities or worst-case behavior.
- Subspace embedding quality for every input matrix.
- Downstream low-rank approximation error.
- Computational runtime.

These require separate analysis or experiments.

---

## 5. Consolidated Findings

| Question | Prediction | Measurement | Conclusion |
|---|---|---|---|
| Does $\zeta$ reduce global Frobenius distortion? | No second-moment dependence on $\zeta$ | 2.039 measured vs. 2.086 predicted; flat over $\zeta=2,\ldots,32$ | Supported for the tested metric |
| Is there an accuracy penalty? | Potentially at very low sparsity | $\zeta=1$ was 8–12 SE worse; $\zeta\geq2$ within 1.3 SE | Penalty detected only at $\zeta=1$ |
| Does rotation change the comparison? | Possible because sparse sign is not rotation-invariant | No effect detected at $n=200$ | Null result within tested design |
| Is sparse sign computationally cheaper? | Expected when sparsity is exploited | Not measured; sketch returned dense | Open |

## Next Steps

1. Preserve sparse matrix structure and benchmark actual construction, multiplication, and memory costs.
2. Investigate whether $\zeta$ affects spectral distortion or tail behavior despite the Frobenius second-moment independence.
3. Revisit orientation sensitivity using a design targeted at detecting alignment-dependent effects.
4. Compare these findings with the sparse embedding and leverage-score literature in Pearce and Martinsson's RNLA reading list.

The main result so far is that **additional nonzeros do not improve the expected squared global Gram distortion under this sketch construction**, while the experiments identify $\zeta=1$ as a distinct low-accuracy regime. The computational benefit of sparse sketches remains a separate, unmeasured question.