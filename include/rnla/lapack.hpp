#pragma once

extern "C" {

// C := alpha * op(A) * op(B) + beta * C, op(X) = X or X^T.
void dgemm_(
    const char* transa,
    const char* transb,
    const int* m,
    const int* n,
    const int* k,
    const double* alpha,
    const double* a,
    const int* lda,
    const double* b,
    const int* ldb,
    const double* beta,
    double* c,
    const int* ldc
);

void dgesdd_(
    const char* jobz,
    const int* m,
    const int* n,
    double* a,
    const int* lda,
    double* s,
    double* u,
    const int* ldu,
    double* vt,
    const int* ldvt,
    double* work,
    const int* lwork,
    int* iwork,
    int* info
);

// QR factorization: R in the upper triangle of A, Householder reflectors below,
// with their scalar factors in tau.
void dgeqrf_(
    const int* m,
    const int* n,
    double* a,
    const int* lda,
    double* tau,
    double* work,
    const int* lwork,
    int* info
);

// Expands the reflectors from dgeqrf_ into an explicit orthonormal Q.
void dorgqr_(
    const int* m,
    const int* n,
    const int* k,
    double* a,
    const int* lda,
    const double* tau,
    double* work,
    const int* lwork,
    int* info
);

}