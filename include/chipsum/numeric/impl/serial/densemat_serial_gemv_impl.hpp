#ifndef __CHIPSUM_DENSEMAT_SERIAL_GEMV_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_SERIAL_GEMV_IMPL_HPP__

#include <vector>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType,typename DT>
// gemv
CHIPSUM_FUNCTION_INLINE void
gemv(const DT &A,
     const ::std::vector<ValueType> &x,
     ::std::vector<ValueType> &b) {

    ::std::size_t M = A.nrow;
    ::std::size_t N = A.ncol;
    for (::std::size_t i = 0; i < M; ++i)
        b[i] = 0;

    for (::std::size_t i = 0; i < M; ++i) {
        for (::std::size_t j = 0; j < N; ++j) {
            b[i] += A.data[i * N + j] * x[j];
        }
    }
}

}
}
}
}
#endif // DENSEMAT_SERIAL_GEMM_IMPL_HPP
