#include "densemat_blas_impl.hpp"

template <>
void ChipSum::Numeric::DenseMat::Mult<double>(blasint M, blasint N, blasint K, double *A, double *B, double *C)
{
    cblas_dgemm(CBLAS_ORDER::CblasRowMajor,
                CBLAS_TRANSPOSE::CblasNoTrans,
                CBLAS_TRANSPOSE::CblasNoTrans,
                M,N,K,1.0,A,K,B,N,0.0,C,N
                );
}
