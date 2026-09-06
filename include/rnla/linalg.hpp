#pragma once

#include <vector>

#include "rnla/matrix.hpp"

namespace rnla {

// Rank-k truncated SVD: A ≈ U * diag(s) * Vt.
struct TruncatedSVD {
  Matrix U;               // m x k
  std::vector<double> s;  // k singular values, descending
  Matrix Vt;              // k x n
};

TruncatedSVD truncated_svd(const Matrix& A, int k);

// All singular values of A, in descending order.
std::vector<double> singular_values(const Matrix& A);

// Orthonormal basis for the column space of A, via Householder QR.
// Assumes A has full column rank.
Matrix orth(const Matrix& A);

// Orthonormal basis for the column space of A via CholeskyQR:
//   G = A^T A,  G = R^T R,  Q = A R^{-1}
// Three BLAS-3 operations, far more parallel than Householder QR — but it
// forms the Gram matrix, so the Cholesky sees cond(A)^2 and the method fails
// once cond(A) exceeds about eps^{-1/2} ~ 1e8.
// Throws if the Gram matrix is not numerically positive definite.
Matrix cholesky_qr(const Matrix& A);

// CholeskyQR applied twice. The first pass leaves Q1 with
// ||Q1^T Q1 - I|| ~ eps*cond(A)^2, which makes cond(Q1) close to 1 whenever
// that quantity is below 1 — so the second pass is well-conditioned and
// returns machine-precision orthogonality.
Matrix cholesky_qr2(const Matrix& A);

// CholeskyQR with a diagonal shift on the Gram matrix: G = A^T A + s*I.
// The shift keeps G numerically positive definite past the point where plain
// CholeskyQR fails, at the cost of a Q that is not orthonormal.
Matrix cholesky_qr_shifted(const Matrix& A, double s);

// Spectral norm (largest singular value). Computes a full SVD, so this is for
// diagnostics and test harnesses, not hot paths.
double norm_2(const Matrix& A);

// Shifted CholeskyQR3 (Fukaya et al. 2020). One shifted pass to condition the
// problem, then two plain passes to restore orthonormality.
Matrix cholesky_qr3(const Matrix& A);

// ||Q^T Q - I||_F. Measures how far Q is from having orthonormal columns.
double orthogonality_error(const Matrix& Q);

// Frobenius norm.
double norm_fro(const Matrix& A);

// ||A - U diag(s) Vt||_F  — how well the factorization reproduces A.
double reconstruction_error(const Matrix& A, const TruncatedSVD& svd);

// Eckart-Young optimum in Frobenius norm: sqrt(sum_{i>k} sigma_i^2).
// The best any rank-k matrix can do, given the true spectrum.
double eckart_young_fro(const std::vector<double>& sigma, int k);

}  // namespace rnla