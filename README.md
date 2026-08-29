# Randomized Numerical Linear Algebra

**Randomized low-rank approximation in C++20: randomized SVD and its variants, measured for accuracy against the Eckart–Young optimum.**

## Problem

Many matrices in scientific computing and machine learning are large but approximately low rank: a covariance matrix, a kernel Gram matrix, a neural network weight matrix. The best rank-$k$ approximation is given by the truncated SVD, but computing it costs $O(mn\min(m,n))$ - you pay for the entire spectrum in order to discard most of it.

Randomized methods invert this. Probe the matrix with $\ell = k + p$ random vectors, orthonormalize the result, and project:

$$Y = A\Omega, \qquad Q = \operatorname{orth}(Y), \qquad B = Q^\top A, \qquad A \approx Q B_k$$

The cost scales with the rank you want, not the rank the matrix has. The obvious question is what accuracy that buys, and the usual answer is a worst-case bound.

We ask a different question: **what actually governs the error in practice, and is it predictable from the spectrum before running anything?**

## Hypothesis

The natural guess is the local gap at the truncation point, $\sigma_{k+1}/\sigma_k$ - a sharply separated top-$k$ subspace should be easy to find, a poorly separated one hard.

**This is measurably false.** Two spectra with gaps of 0.9048 and 0.9070 give ratios to optimal of 1.199 and 1.319, and the error curve is *non-monotonic* in the gap: it rises, peaks near 0.93, then falls back to exactly 1.000 at a perfectly flat spectrum.

The better candidate is the **effective tail dimension**

$$\mathrm{sr}_{\text{tail}} = \frac{\sum_{j>k}\sigma_j^2}{\sigma_{k+1}^2},$$

which counts how many directions the tail effectively occupies. Both extremes are easy - a tail of dimension 1 means the top-$k$ subspace is sharply defined, and a very large tail dimension means every $k$-dimensional subspace is equally good - with maximum difficulty in between, where the tail dimension is comparable to the sketch width $\ell$.

A second finding follows, and is arguably the more useful one: **power iteration changes which variable governs the error.** At $q = 0$ the curve is non-monotonic in the gap; at $q \ge 1$ it is monotone increasing, with both spectral families interleaving in gap order. The non-monotonic structure is a $q = 0$ phenomenon.

## Methodology

Every result is a **ratio to the Eckart–Young optimum**, not an absolute error. Test matrices are built as $A = U\,\mathrm{diag}(\sigma)\,V^\top$ with Haar-random orthonormal factors, so the spectrum is exactly as prescribed (verified to $2\times10^{-15}$) and the optimum

$$\|A - A_k\|_F = \Big(\sum_{i>k}\sigma_i^2\Big)^{1/2}$$

is available in closed form. Nothing is measured against another implementation.

The comparison set:

- Randomized SVD with oversampling $p$ (HMT Alg. 4.1 / 5.1)
- Power iteration with re-orthonormalization, $q$ passes (HMT Alg. 4.4)
- Randomized block Krylov (Musco & Musco) - *next*
- Deterministic truncated SVD via LAPACK `dgesdd` as ground truth

Spectra are polynomial ($\sigma_i = (i+1)^{-\alpha}$), exponential ($\sigma_i = e^{-\alpha i}$), and flat, swept across $\alpha$ so that the gap and the tail dimension take a range of values across both families.

Two results worth singling out:

**Flat spectra are exactly optimal, provably.** The projection loses $r - \ell$ and truncating $B$ from rank $\ell$ to $k$ loses $\ell - k$; the terms are orthogonal and sum to exactly $r - k$. Ratio $= 1$ for every draw of $\Omega$, confirmed to nine decimals with zero variance across ten seeds.

**Naive power iteration fails at a predictable point.** The dynamic range across the sketch grows as $(\sigma_1/\sigma_\ell)^{2q+1}$; once it exceeds $1/\varepsilon$ the small directions fall below the largest one's rounding error. The failure at $q = 3$ was predicted from the spectrum before it was observed, and $q = 2$ was already silently 2.5× degraded.

Full derivations in [`derivations/`](derivations/).

## Implementation

- **C++20 throughout.** LAPACK/BLAS underneath - the contribution is algorithmic and empirical, not a BLAS reimplementation.
- **Fortran symbols declared by hand** in `include/rnla/lapack.hpp` rather than pulled from a vendor header. Accelerate, OpenBLAS and MKL share the Fortran ABI but ship conflicting C headers; declaring the dozen symbols directly meant the Linux port required zero source changes.
- **Column-major storage**, because LAPACK requires it. Storing row-major and transposing at each call costs a full copy per factorization, which would silently dominate any timing.
- **Every random draw is explicitly seeded.** No defaults, no `random_device`. Every result is a sample from a distribution and must be reproducible.
- **CUDA (planned)** - only once profiling identifies a bottleneck.

## Building

```bash
cmake --preset default
cmake --build build
ctest --test-dir build --output-on-failure
./build/bench_spectrum
```

Requires CMake ≥ 3.21 and a C++20 compiler. BLAS/LAPACK is Apple Accelerate on macOS and the system LAPACK on Linux (`liblapack-dev libblas-dev`); CMake selects automatically and reports which backend it found.

An `asan` preset builds with AddressSanitizer and UndefinedBehaviorSanitizer for running the test suite.

## Tech stack

```
C++20                 no exceptions to the standard
CMake + presets       default (Release) and asan (Debug + sanitizers)
vcpkg                 manifest mode, per-project dependencies
LAPACK / BLAS         Accelerate on macOS, system LAPACK on Linux
CTest                 hand-rolled check()/close(), no framework yet
clang-format          Google style
```

### Compute

- **NVIDIA RTX 5060 Ti**, 16 GB, via WSL2 (CUDA 13.0)
- **CARC** (carc.unm.edu) for anything at scale
- MacBook M2 for development; Accelerate is fast on CPU, **no CUDA**

## Repository structure

```
randomized-numerical-linear-algebra/
├── README.md
├── CMakeLists.txt
├── CMakePresets.json               # default + asan
├── vcpkg.json
├── include/rnla/
│   ├── matrix.hpp                  # column-major dense matrix
│   ├── lapack.hpp                  # hand-written extern "C" declarations
│   ├── linalg.hpp                  # SVD, QR, norms, error metrics
│   ├── random.hpp                  # seeded Gaussian sketch
│   ├── test_matrices.hpp           # prescribed-spectrum generator
│   └── rsvd.hpp                    # randomized SVD + variants
├── src/                            # implementations, one per header
├── tests/
│   └── test_linalg.cpp             # exact identities + algorithmic invariants
├── benchmarks/
│   └── bench_spectrum.cpp          # reports, does not assert
├── derivations/                    # the mathematics, from first principles
│   └── notes/                      # standalone results from experiments
└── docs/
    └── project-map.md              # status, decisions, open directions
```

**Test suite:**
- Deterministic SVD attains the Eckart–Young optimum to machine precision.
- Generated matrices have exactly the requested singular values.
- Randomized methods never beat Eckart–Young - a theorem, so a strong bug detector.
- Power iteration error is monotone decreasing in $q$, guarding the re-orthonormalization.

Benchmarks are deliberately not CTest targets. Tests assert; benchmarks report.

## Status and roadmap

V1 complete as of Aug 2026: randomized SVD with oversampling and stable power iteration, characterized across 11 spectra × 3 values of $q$ × 10 seeds. Full plan, hypotheses, and known debts in [`docs/project-map.md`](docs/project-map.md).

| Phase | Content | Target deliverable | Date |
|---|---|---|---|
| **0. Foundation** | Matrix type, LAPACK bindings, two-platform build. | Clean build on macOS and Linux. | Aug 2026 |
| **1. Baseline and rSVD** | Truncated SVD, `orth`, Gaussian sketch, spectrum generator, rSVD with oversampling and power iteration. | **H1–H5 tested.** Spectrum sweep. | Aug 2026 |
| **2. Block Krylov** | Musco–Musco; comparison at matched passes over $A$. | **H6 tested.** | September |
| **3. Orthogonalization** | CholeskyQR / CholeskyQR2 / shifted CholeskyQR3 / TSQR. | **H7 tested.** | October |
| **4. Sketch operators** | SRHT, sparse sign; first timing harness. | **H8 tested.** | November |
| **5. Performance** | Profiling, roofline, GPU where the profiler points. | Runtime breakdown. | Winter break |
| **6. Write-up** | Consolidate derivations. | Readable end to end. | Spring 2027 |

**Note. No performance claims until a fair deterministic baseline exists** - the current `truncated_svd` computes the full SVD and discards the tail, which flatters every randomized method compared against it.

## References

- **Halko, Martinsson & Tropp (2011)**, *Finding Structure with Randomness*, SIAM Review.
- **Martinsson & Tropp (2020)**, *Randomized Numerical Linear Algebra: Foundations and Algorithms*, Acta Numerica.
- **Musco & Musco (2015)**, *Randomized Block Krylov Methods for Stronger and Faster Approximate SVD*, NIPS.
- **Tropp, Yurtsever, Udell & Cevher (2017)**, *Practical Sketching Algorithms for Low-Rank Matrix Approximation*, SIMAX.
- **Clarkson & Woodruff (2013)**, *Low-Rank Approximation and Regression in Input Sparsity Time*, STOC.
- **Frangella, Tropp & Udell (2023)**, *Randomized Nyström Preconditioning*, SIMAX.
- **Murray et al. (2023)**, *Randomized Numerical Linear Algebra: A Perspective on the Field with an Eye to Software*.
- **Fukaya et al. (2020)**, *Shifted CholeskyQR3 for Computing the QR Factorization of Ill-Conditioned Matrices*, SISC.
- **Golub & Van Loan**, *Matrix Computations*, chs. 5, 8, 10.
- **Higham**, *Accuracy and Stability of Numerical Algorithms*, chs. 19–20.