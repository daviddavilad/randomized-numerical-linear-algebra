# RNLA Project Map - Randomized Low-Rank Approximation

**Status:** V2 complete. Block Krylov measured; performance work not started.
**Author:** David Dávila
**Date:** 1 Sep 2026

---

## 1. Title

### Working titles (tentative)

1. *What Governs the Accuracy of Randomized Low-Rank Approximation?* - most acurate, currently the whole project.
2. *Spectral Geometry and the Cost of Randomized Sketching* - use if the tail-dimension result firms up into the central claim.
3. *Gap-Independence in Practice: Block Krylov vs. Power Iteration* - use if the V2 comparison turns out to be the strongest finding.

Start with (1) and reevaluate once V2 lands.

### The research question

Randomized SVD comes with worst-case error bounds. In practice it usually does much better than the bound, sometimes exactly optimally, and occasionally fails outright - for numerical rather than statistical reasons.

This project characterizes the actual behaviour: which spectral features predict the error, when each algorithmic variant is worth its cost, and where the floating-point limits bite. Every claim is a ratio to the Eckart–Young optimum,

$$\min_{\mathrm{rank}(B)\le k}\|A-B\|_F = \Big(\sum_{i>k}\sigma_i^2\Big)^{1/2},$$

measured on matrices whose spectra are known exactly by construction rather than estimated. Nothing here is measured against another implementation.

The output is a characterization rather than a single thesis: each result derived where a derivation exists, measured where it doesn't, and predicted in advance where possible.

### Project complementarity

This repository supplies infrastructure and methods to [scalable-kernel-learning](https://github.com/daviddavilad/scalable-kernel-learning) - the prescribed-spectrum generator and the randomized range finder - but asks a different question.

> **RNLA:** how well does a sketch reproduce a *matrix*?
> **SKL:** how well does an approximation reproduce *predictions*?
> **parallel-cg-ridge:** how do we scale the exact solve, without approximating at all?

The distinction between the first two is the point rather than a technicality. SKL's central claim is that matrix error is the wrong target for a learning problem: two approximations with equal spectral error can behave differently downstream. RNLA builds the instrument; SKL asks whether the instrument measures the right thing.

---

## 2. Hypotheses and Findings

| # | Hypothesis | Status | Falsification |
|---|---|---|---|
| **H1** | Error is governed by the local gap $\sigma_{k+1}/\sigma_k$. | **Rejected** | Two spectra with matched gaps (0.9048, 0.9070) gave 1.199 and 1.319. |
| **H2** | Error is governed by effective tail dimension $`\mathrm{sr}_\mathrm{tail} = (\sum_{j>k}\sigma_j^2)/\sigma_{k+1}^2`$; difficulty peaks where it is comparable to the sketch width $`\ell`$. | **Supported, not decisive** | A spectrum pair with matched $`\mathrm{sr}_\mathrm{tail}`$ and different error. |
| **H3** | Flat spectra give ratio exactly 1 for every draw of $\Omega$. | **Proved and confirmed** | Any seed-to-seed variance, or ratio $\ne 1$. |
| **H4** | Naive power iteration fails once $(\sigma_1/\sigma_\ell)^{2q+1}$ exceeds $1/\varepsilon$. | **Confirmed, predicted in advance** | Failure at a $q$ unrelated to that threshold. |
| **H5** | Power iteration restores gap-monotonicity: at $q \ge 1$ the error curve is monotone in the gap. | **Observed at $q = 1, 2$** | An inversion at $q \ge 1$, or non-monotonicity returning at larger $q$. |
| **H6** | Block Krylov's gap-independent bound shows up empirically as a *flatter* error curve across gaps. | **Rejected** | The spread across polynomial rows widened from 104× to 1708×; the improvement shrinks as the gap approaches 1. The prediction mistranslated the theorem - see `notes/block-krylov-prediction.md`. |
| **H7** | CholeskyQR2 recovers Householder-quality orthogonality at substantially lower cost, up to a predictable conditioning threshold. | **Open - V3** | CholeskyQR2 fails where CholeskyQR does, or shows no speed advantage. |
| **H8** | Sketch operator choice (Gaussian / SRHT / sparse sign) changes cost but not accuracy, outside pathological spectra. | **Open - V4** | Measurable accuracy differences at matched $\ell$. |

H2 is the central open claim. H7 is the next test.

---

## 3. Scope

### In scope

**Baselines (fp64, ground truth)**
- Deterministic truncated SVD via LAPACK `dgesdd`
- Eckart–Young optimum computed in closed form from the prescribed spectrum

**Algorithms**
- Randomized SVD with oversampling (Halko–Martinsson–Tropp Alg. 4.1 / 5.1)
- Power iteration with re-orthonormalization (HMT Alg. 4.4)
- Randomized block Krylov (Musco & Musco 2015)
- Orthogonalization: Householder QR; CholeskyQR, CholeskyQR2, shifted CholeskyQR3, TSQR - *next*
- Sketch operators: Gaussian; SRHT, sparse sign / CountSketch - *planned*

**Diagnostics**
- Ratio to the Eckart–Young optimum in Frobenius norm
- $\|Q^\top Q - I\|_F$ as an orthogonality measure
- Spectral gap, effective tail dimension, dynamic range across the sketch

### Potential extensions

| Extension | Explanation |
|---|---|
| **Randomized Nyström preconditioning** | Frangella, Tropp & Udell (2023): use a randomized low-rank approximation to deflate the top eigenvalues of an SPD system, collapsing CG's iteration count. Directly connects this repo to the MATH 471 solver, and to Prof. Schroder's work on randomization in multigrid. The strongest single extension. |
| **Mixed precision** | Sketch in fp32/bf16 and refine in fp64. Interacts with H7: CholeskyQR in low precision is far more fragile. Higham & Mary, Carson & Higham. |
| **Adaptive rank selection** | HMT §4.4 gives a posterior error estimator, so `target_error → k` rather than a fixed $k$. Reproduction rather than contribution, but genuinely useful. |
| **Matrix-free operators** | $A$ available only as a matvec. Required for anything at scale, and the doorway to Krylov methods and PDE operators. |
| **GPU / CUDA** | Only once profiling identifies a bottleneck. The sketch phase is memory-bandwidth-bound and the QR is compute-bound - that asymmetry should drive the work. |
| **ML-relevant matrices** | Neural network weight matrices, attention matrices, Hessian approximations, LoRA-style updates. Tests whether spectral error is even the right objective for downstream ML behaviour. |

### Out of scope

- **Rewriting GEMM or QR from scratch.** Losing to OpenBLAS is not a result. The contribution is algorithmic and empirical, not a BLAS reimplementation.
- **Competing with cuSOLVER on raw throughput.** Also not a result.
- **Reimplementing RandLAPACK/RandBLAS.** A funded standardization effort already exists (Murray et al. 2023); position beside it, not against it.
- **Sparse matrices.** Dense only, for now. A different set of bottlenecks.
- **Non-symmetric eigenproblems.** Different project.

---

## 4. Experimental design

### Test matrices

$A = U\,\mathrm{diag}(\sigma)\,V^\top$ with Haar-random orthonormal $U, V$, so the spectrum is exactly as prescribed (verified to $2\times10^{-15}$).

| Spectrum | Form | Purpose |
|---|---|---|
| Polynomial | $\sigma_i = (i+1)^{-\alpha}$ | Heavy tail; the hard regime for randomized methods |
| Exponential | $\sigma_i = e^{-\alpha i}$ | Fast decay; easy regime |
| Flat | $\sigma_i = 1$ | Degenerate case with an exact closed-form answer |

Current sweep: 11 spectra × 10 seeds at $`m = 300`$, $`n = 200`$, $`k = 20`$, $`p = 10`$ - power iteration at $`q = 0, 1, 2`$ and block Krylov at $`q = 1`$.

### Metrics

*Accuracy:* ratio to the Eckart–Young optimum (mean and max over seeds); $\|Q^\top Q - I\|_F$.

*Computational (not yet instrumented):* wall-clock separated by phase (sketch, orthogonalization, projection, small SVD); peak memory; passes over $A$; analytic FLOP count.

*Diagnostic:* $`\sigma_{k+1}/\sigma_k`$; $`\mathrm{sr}_\mathrm{tail}`$; dynamic range $`(\sigma_1/\sigma_\ell)^{2q+1}`$.

### Important research figures (deliverables)

1. Ratio-to-optimal vs. gap, one line per $q$ - **the money plot.** If H5 holds, $q\ge1$ curves are monotone where $q=0$ is not.
2. Ratio-to-optimal vs. $\mathrm{sr}_{\text{tail}}$, families overlaid - the H2 plot. If H2 holds, the two families collapse onto one curve.
3. Block Krylov vs. power iteration at matched passes over $A$ - the H6 plot.
4. $\|Q^\top Q - I\|_F$ vs. $\kappa$, one line per orthogonalization scheme, with the CholeskyQR breakdown marked - the H7 plot.
5. Phase-wise runtime breakdown vs. $(m, n, k)$ - where the time actually goes.
6. Accuracy per unit cost, all variants overlaid - the Pareto frontier.

---

## 5. Tech stack

### C++ (the whole project)

```
C++20                 no exceptions to the standard
CMake + presets       default (Release) and asan (Debug + sanitizers)
vcpkg                 manifest mode, per-project dependencies
LAPACK / BLAS         Accelerate on macOS, system LAPACK on Linux
CTest                 hand-rolled check()/close(), no framework yet
clang-format          Google style
```

Fortran symbols are declared by hand in `lapack.hpp` rather than pulled from a vendor header, because Accelerate/OpenBLAS/MKL share the Fortran ABI but ship conflicting C headers. This is why the Linux port required zero source changes.

### CUDA (planned)

Only after CPU profiling identifies a bottleneck. Custom kernels benchmarked against cuBLAS/cuSOLVER formulations, never as a replacement for them.

### Compute

- **NVIDIA RTX 5060 Ti**, 16 GB - via WSL2, CUDA 13.0
- **CARC** (carc.unm.edu) - for anything needing datacenter hardware or scale
- **MacBook M2** for development; Accelerate uses the AMX coprocessor and is fast for CPU work, but **no CUDA**

### Writing

`derivations/` in GitHub-rendered Markdown so it reads in the browser. LaTeX reserved for a `paper/` directory if one is ever warranted.

---

## 6. Repository structure

```
randomized-numerical-linear-algebra/
├── README.md                       # problem, hypothesis, method, build
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
│   └── bench_spectrum.cpp          # deliberately NOT a CTest target
├── derivations/
│   ├── 01-svd-and-low-rank.md      # the linear development
│   ├── 02-randomized-range-finding.md
│   ├── 03-randomized-svd.md
│   └── notes/                      # standalone results from experiments
│       ├── flat-spectrum-exact.md
│       ├── power-iteration-stability.md
│       ├── power-iteration-gap-lookup.md
│       └── block-krylov-prediction.md
└── docs/
    └── project-map.md
```

**Test suite:**
- Deterministic SVD attains the Eckart–Young optimum to machine precision.
- Generated matrices have exactly the requested singular values.
- Randomized methods never beat Eckart–Young (a theorem, so a strong bug detector).
- Power iteration error is monotone decreasing in $q$ (guards the re-orthonormalization).

---

## 7. Roadmap

| Phase | Content | Target deliverable | Date |
|---|---|---|---|
| **0. Foundation** | Matrix type, LAPACK bindings, build on two platforms. | Clean build on macOS and Linux; `dgemm` verified by hand. | Aug 2026 |
| **1. Baseline and rSVD** | Truncated SVD, `orth`, Gaussian sketch, prescribed-spectrum generator, rSVD with oversampling and stable power iteration. | **H1–H5 tested.** Spectrum sweep across 11 spectra × 3 $q$ × 10 seeds. | Aug 2026 |
| **2. Block Krylov** | Musco–Musco block Krylov; comparison at matched passes over $A$. | **H6 tested.** Figure 3. | September |
| **3. Orthogonalization** | CholeskyQR, CholeskyQR2, shifted CholeskyQR3, TSQR; conditioning study. | **H7 tested.** Figure 4. | October |
| **4. Sketch operators** | SRHT and sparse sign against Gaussian; first timing harness. | **H8 tested.** Fair deterministic baseline exists. | November |
| **5. Performance** | Phase-wise profiling, roofline analysis, then GPU only where the profiler points. | Figures 5 and 6. | Winter break |
| **6. Write-up** | Consolidate derivations; decide whether a paper is warranted. | Repo readable end to end by a stranger. | Spring 2027 |

**Note. No performance claims until a fair deterministic baseline exists.** The current `truncated_svd` computes the full SVD and discards the tail, which flatters every randomized method it is compared against.

### Known debts

- No timing instrumentation of any kind; every result so far is accuracy only.
- `truncated_svd` is not a fair timing baseline (above).
- Measurement floor at $\sim10^{-8}$ excess: the optimum is analytic, the achieved error numerical, and below that scale they cross. Ratios under 1 appear and are noise.
- Cross-platform agreement verified only to printed precision (4–5 decimals).
- `orth` requires $m \ge n$; wide inputs unhandled.

---

## 8. Readings and resources

### Base understanding

- **Halko, Martinsson & Tropp (2011)**, *Finding Structure with Randomness*, SIAM Review. **The** reference; §4–5 are the algorithms implemented here, §10 the oversampling analysis.
- **Martinsson & Tropp (2020)**, *Randomized Numerical Linear Algebra: Foundations and Algorithms*, Acta Numerica. The modern survey.
- **Eckart & Young (1936)** / **Mirsky (1960)** - the optimality theorem every measurement is relative to.

### Core readings

- **Musco & Musco (2015)**, *Randomized Block Krylov Methods for Stronger and Faster Approximate SVD*, NIPS. **Next implementation.**
- **Tropp, Yurtsever, Udell & Cevher (2017)**, *Practical Sketching Algorithms for Low-Rank Matrix Approximation*, SIMAX. Single-pass methods.
- **Clarkson & Woodruff (2013)**, *Low-Rank Approximation and Regression in Input Sparsity Time*, STOC. Sparse sketching.
- **Frangella, Tropp & Udell (2023)**, *Randomized Nyström Preconditioning*, SIMAX. **The bridge to MATH 471 and to Prof. Schroder's interests.**
- **Murray et al. (2023)**, *Randomized Numerical Linear Algebra: A Perspective on the Field with an Eye to Software*. Know what already exists.

### Background

- **Golub & Van Loan**, *Matrix Computations*, chs. 5, 8, 10.
- **Higham**, *Accuracy and Stability of Numerical Algorithms*, chs. 19–20.
- **Fukaya et al. (2020)**, *Shifted CholeskyQR3*. For the orthogonalization phase.
- **Demmel et al. (2012)**, *Communication-Avoiding Parallel and Sequential QR*. TSQR; also relevant to MATH 471.

---

## 9. Open questions

This project is self-directed, so this section holds questions to resolve rather than questions for an advisor.

1. **Is $\mathrm{sr}_{\text{tail}}$ actually the right variable, or is it a proxy?** H2 is supported but not decisive. A spectrum pair with matched $\mathrm{sr}_{\text{tail}}$ and different error would settle it either way.
2. **Where does block Krylov stop being worth its memory?** $(q{+}1)\ell$ columns rather than $\ell$. At what $q$ and what spectrum does the accuracy gain stop paying for the storage?
3. **Does the $q=0$ non-monotonicity have a closed form?** The flat case was derivable exactly. Whether the peak location is predictable from the spectrum is unknown.
4. **Worth raising with Prof. Schroder:** randomized Nyström preconditioning sits at the intersection of this repo, the MATH 471 solver, and his own work on randomization and multigrid. Best entry point for a supervised extension.
5. **What would actually test Musco–Musco?** Their gap-independence concerns the iteration count needed to reach a target $`\varepsilon`$, not the error at fixed $`q`$. Invert the benchmark: fix a target ratio (say 1.001), and for each spectrum report the smallest $`q`$ at which each method reaches it. Gap-independence predicts block Krylov's required $`q`$ is roughly constant across spectra while power iteration's grows as the gap approaches 1.

---

## 10. Next task

Implement CholeskyQR and CholeskyQR2 as alternatives to Householder QR in `orth`, and measure $`\|Q^\top Q - I\|_F`$ against $`\kappa(Y)`$ for each.

CholeskyQR is three BLAS-3 operations against Householder's sequential reflector sweep, but it forms $`Y^\top Y`$ and so squares the condition number. Prediction to record before running: plain CholeskyQR loses orthogonality once $`\kappa(Y) > \varepsilon^{-1/2} \approx 10^8`$, and CholeskyQR2 recovers machine-precision orthogonality up to roughly that same threshold.