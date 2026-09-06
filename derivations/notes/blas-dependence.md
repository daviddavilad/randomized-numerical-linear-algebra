# A timing result is a property of the BLAS

*2026-09-05. Found while measuring the cost half of H7.*

## Results

The orthogonalization timing benchmark was first run against Ubuntu's default BLAS, which is reference netlib - a naive triple-loop implementation with no blocking, no vectorization, no cache tuning. Installing `libopenblas-dev` and re-running **inverted the conclusion**.

| method | reference BLAS | OpenBLAS (1 thread) | speedup |
|---|---|---|---|
| Householder | 130 µs | 41.5 µs | 3.1× |
| CholeskyQR | 117 µs | 20.3 µs | 5.8× |
| CholeskyQR2 | 235 µs | 37.1 µs | 6.3× |

Under reference BLAS, CholeskyQR2 is **1.8× slower** than Householder. Under OpenBLAS it is **1.12× faster** under the same conditions (same code, same machine, same inputs).

## Explanation

The CholeskyQR routines gained roughly twice as much from optimized BLAS as Householder did - the arithmetic-intensity argument appearing as a number rather than as prose.

CholeskyQR is three BLAS-3 calls (`dgemm`, `dpotrf`, `dtrsm`), the routines implementers optimize hardest and the ones with enough arithmetic per byte moved to reward blocking. Householder QR is a sequence of rank-1 updates with a dependency chain between them; the blocked variant recovers some of this but not all.

**The method's advantage is entirely an advantage over a good BLAS.** Measured against an unoptimized one it does not exist. Under OpenBLAS the times track the flop counts (CholeskyQR ~0.4 MFLOP vs Householder ~0.7, predicting 1.75×, measured 2.05×). Under reference BLAS they do not. Running near the flop count is what "optimized" means.

## Which results are BLAS-dependent

| quantity | reference → OpenBLAS |
|---|---|
| rSVD/optimal ratio | identical to all 12 printed digits |
| $`\|Q^\top Q - I\|_F`$, Householder | 2.6e-15 → 1.4e-15 |
| $`\|Q^\top Q - I\|_F`$, CholeskyQR at $`\kappa=10^2`$ | 6.5e-13 → 2.0e-13 |
| constant $`c`$ in $`c\,\varepsilon\kappa^2`$ | ~0.1 → ~0.03 |

$`O(1)`$ quantities are unaffected: the ratio to the Eckart–Young optimum is a real number that rounding perturbs at the 16th digit, far below anything reported. $`O(\varepsilon)`$ stability diagnostics *are* made of rounding, so accumulation order determines them. The $`\varepsilon\kappa^2`$ scaling law survives; only its constant moves.

## Consequences

1. Every timing number in this repository is OpenBLAS 0.3.26, single-threaded, on one machine. Accelerate and MKL will differ.
2. Threading is a net loss at this size - 54 µs vs 41.5 µs for Householder. Coordination costs more than the parallelism buys at $`m=200`$, $`n=30`$. Benchmarks run with `OPENBLAS_NUM_THREADS=1`.
3. Accuracy results need no such caveat.