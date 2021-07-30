#ifndef __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__

#include <type_traits>
#include <cblas.h>


#include "../numeric_traits.hpp"
#include "../../chipsum_macro.h"


namespace ChipSum {

namespace  Numeric {


template<typename ...Props>
struct DenseMatrix_Traits<double,blasint,ChipSum::Backend::OpenBlas,Props...>
        : public Operator_Traits<double,blasint,ChipSum::Backend::OpenBlas,Props...>{



    using matrix_type = double*;

    using size_type = blasint;


};



namespace  DenseMat
{


CHIPSUM_FUNCTION_INLINE void Mult(blasint M,blasint N,blasint K,double* A,double* B,double* C)
{
    cblas_dgemm(CBLAS_ORDER::CblasRowMajor,
                CBLAS_TRANSPOSE::CblasNoTrans,
                CBLAS_TRANSPOSE::CblasNoTrans,
                M,N,K,1.0,A,K,B,N,0.0,C,N
                );
}



} // namespace DenseMat
} // namespace Numeric
} // namespace ChipSum
#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
