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

}  // namespace rnla