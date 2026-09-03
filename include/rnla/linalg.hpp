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