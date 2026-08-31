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

void set_block(Matrix& dst, const Matrix& src, int col_offset) {
  assert(dst.rows() == src.rows());
  assert(col_offset >= 0 && col_offset + src.cols() <= dst.cols());

  // Copy, j outer over src.cols(), i inner over src.rows()
  for (int j = 0; j < src.cols(); ++j) {
    for (int i = 0; i < src.rows(); ++i) {
      dst(i, col_offset + j) = src(i, j);
    }
  }
}

}  // namespace rnla