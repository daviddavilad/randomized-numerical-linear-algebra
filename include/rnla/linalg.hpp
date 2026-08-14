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

// ||Q^T Q - I||_F. Measures how far Q is from having orthonormal columns.
double orthogonality_error(const Matrix& Q);

// Frobenius norm.
double norm_fro(const Matrix& A);

}  // namespace rnla