#include "rnla/matrix.hpp"

#include <cassert>

#include "rnla/lapack.hpp"

namespace rnla {

Matrix matmul(const Matrix& A, const Matrix& B, bool transA, bool transB) {
  const int m = transA ? A.cols() : A.rows();
  const int n = transB ? B.rows() : B.cols();
  const int k = transA ? A.rows() : A.cols();
  const int kb = transB ? B.cols() : B.rows();
  assert(k == kb && "inner dimensions must agree");
  (void)kb;  // silences an unused-variable warning in release builds

  Matrix C(m, n);
  const char ta = transA ? 'T' : 'N';
  const char tb = transB ? 'T' : 'N';
  const double alpha = 1.0, beta = 0.0;
  const int lda = A.ld(), ldb = B.ld(), ldc = C.ld();

  dgemm_(&ta, &tb, &m, &n, &k, &alpha, A.data(), &lda, B.data(), &ldb, &beta, C.data(), &ldc);
  return C;
}

}  // namespace rnla