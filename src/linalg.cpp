#include "rnla/linalg.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <string>
#include <cmath>

#include "rnla/lapack.hpp"

namespace rnla {

std::vector<double> singular_values(const Matrix& A) {
  Matrix work_copy = A;

  const int m = A.rows();
  const int n = A.cols();
  const int lda = A.ld();
  const int mn = std::min(m, n);

  const char jobz = 'N';
  int ldu = 1, ldvt = 1;

  std::vector<double> s(mn);
  std::vector<int> iwork(8 * mn);
  int info = 0;

  // --- workspace query -------------------------------------------------
  double wkopt = 0.0;
  int lwork = -1;
  // Call dgesdd_ with &wkopt as work and lwork = -1
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), nullptr, &ldu, nullptr, &ldvt, &wkopt, &lwork, iwork.data(), &info);

  lwork = static_cast<int>(wkopt);
  std::vector<double> work(std::max(lwork, 1));

  // --- real call -------------------------------------------------------
  // Call dgesdd_ again, with work.data() and the real lwork
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), nullptr, &ldu, nullptr, &ldvt, work.data(), &lwork, iwork.data(), &info);

  if (info != 0)
    throw std::runtime_error("dgesdd failed, info = " + std::to_string(info));

  return s;
}

TruncatedSVD truncated_svd(const Matrix& A, int k) {
  Matrix work_copy = A;

  const int m = A.rows();
  const int n = A.cols();
  const int lda = A.ld();
  const int mn = std::min(m, n);

  if (k < 0 || k > mn)
    throw std::invalid_argument("k must be in [0, min(m, n)]");

  const char jobz = 'S';

  Matrix Ufull(m, mn);
  Matrix Vtfull(mn, n);
  const int ldu = Ufull.ld();
  const int ldvt = Vtfull.ld();

  std::vector<double> s(mn);
  std::vector<int> iwork(8 * mn);
  int info = 0;

  double wkopt = 0.0;
  int lwork = -1;
  // Workspace query
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), Ufull.data(), &ldu, Vtfull.data(), &ldvt, &wkopt, &lwork, iwork.data(), &info);

  lwork = static_cast<int>(wkopt);
  std::vector<double> work(std::max(lwork, 1));

  // Real call
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), Ufull.data(), &ldu, Vtfull.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

  if (info != 0)
    throw std::runtime_error("dgesdd failed, info = " + std::to_string(info));

  TruncatedSVD out;
  out.U = Matrix(m, k);
  out.Vt = Matrix(k, n);
  out.s.assign(s.begin(), s.begin() + k);

  // Copy first k columns of Ufull into out.U
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < m; ++i) out.U(i, j) = Ufull(i, j);

  // Copy first k rows of Vtfull into out.Vt
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < k; ++i) out.Vt(i, j) = Vtfull(i, j);

  return out;
}

Matrix orth(const Matrix& A) {
  Matrix Q = A;

  const int m = A.rows();
  const int n = A.cols();
  const int lda = Q.ld();
  const int k = std::min(m, n);

  if (m < n)
    throw std::invalid_argument("orth requires m >= n (tall or square)");

  std::vector<double> tau(k);
  int info = 0;

  // --- dgeqrf: factorize ------------------------------------------------
  double wkopt = 0.0;
  int lwork = -1;
  // dgeqrf_ workspace query
  dgeqrf_(&m, &n, Q.data(), &lda, tau.data(), &wkopt, &lwork, &info);

  lwork = static_cast<int>(wkopt);
  std::vector<double> work(std::max(lwork, 1));

  // dgeqrf_ real call
  dgeqrf_(&m, &n, Q.data(), &lda, tau.data(), work.data(), &lwork, &info);

  if (info != 0)
    throw std::runtime_error("dgeqrf failed, info = " + std::to_string(info));

  // --- dorgqr: expand the reflectors into Q -----------------------------
  lwork = -1;
  // dorgqr_ workspace query
  dorgqr_(&m, &n, &k, Q.data(), &lda, tau.data(), &wkopt, &lwork, &info);

  lwork = static_cast<int>(wkopt);
  work.assign(std::max(lwork, 1), 0.0);

  // dorgqr_ real call
  dorgqr_(&m, &n, &k, Q.data(), &lda, tau.data(), work.data(), &lwork, &info);

  if (info != 0)
    throw std::runtime_error("dorgqr failed, info = " + std::to_string(info));

  return Q;
}

Matrix cholesky_qr(const Matrix& A) {
  const int m = A.rows();
  const int n = A.cols();
  if (m < n)
    throw std::invalid_argument("cholesky_qr requires m >= n");
  
  Matrix G = matmul(A, A, true, false);   // A^T A, n x n
  const int ldg = G.ld();

  const char uplo = 'U';
  int info = 0;
  dpotrf_(&uplo, &n, G.data(), &ldg, &info);

  if (info < 0)
    throw std::runtime_error("dpotrf: invalid argument " + std::to_string(-info));
  if (info > 0)
    throw std::runtime_error(
        "cholesky_qr: Gram matrix lost positive definiteness at leading minor " +
        std::to_string(info) + " — cond(A) likely exceeds eps^{-1/2}");

  Matrix Q = A;   // copy: A is const, and dtrsm overwrites its operand

  // Q = Q * R^{-1}, i.e. solve X * R = Q for X
  const char side = 'R';
  const char trans = 'N';
  const char diag = 'N';
  const double alpha = 1.0;
  const int ldq = Q.ld();
  dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, G.data(), &ldg, Q.data(), &ldq);

  return Q;
}

double norm_fro(const Matrix& A) {
  double acc = 0.0;
  // Accumulate squares of every entry, j outer
  for (int j = 0; j < A.cols(); ++j)
    for (int i = 0; i < A.rows(); ++i)
      acc += A(i, j) * A(i, j);
  return std::sqrt(acc);
}

double orthogonality_error(const Matrix& Q) {
  Matrix G = matmul(Q, Q, true, false);   // Q^T Q, k x k
  // Subtract 1 from each diagonal entry
  for (int i = 0; i < G.rows(); ++i)
    G(i, i) -= 1.0;
  return norm_fro(G);
}

double eckart_young_fro(const std::vector<double>& sigma, int k) {
  double tail = 0.0;
  for (std::size_t i = static_cast<std::size_t>(k); i < sigma.size(); ++i)
    tail += sigma[i] * sigma[i];
  return std::sqrt(tail);
}

double reconstruction_error(const Matrix& A, const TruncatedSVD& svd) {
  const int k = static_cast<int>(svd.s.size());

  Matrix Us = svd.U;
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < Us.rows(); ++i) Us(i, j) *= svd.s[j];

  Matrix Ak = matmul(Us, svd.Vt);

  Matrix R(A.rows(), A.cols());
  for (int j = 0; j < A.cols(); ++j)
    for (int i = 0; i < A.rows(); ++i) R(i, j) = A(i, j) - Ak(i, j);

  return norm_fro(R);
}

}  // namespace rnla