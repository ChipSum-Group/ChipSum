#ifndef __CHIPSUM_DENSEMAT_SERIAL_GEMM_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_SERIAL_GEMM_IMPL_HPP__

#include <vector>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename DT>
// gemm
CHIPSUM_FUNCTION_INLINE void
gemm(const DT &A,
     const DT &B,
     DT &C) {


    ::std::size_t M = A.nrow;
    ::std::size_t N = B.ncol;
    ::std::size_t K = A.ncol;

    C.nrow = M;
    C.ncol = N;


    for (::std::size_t i = 0; i < C.data.size(); ++i)
        C.data[i] = 0;

    for (::std::size_t m = 0; m < M; ++m) {
        for (::std::size_t k = 0;  k< K; ++k) {
            auto Aik = A.data[m * K + k];
            for (::std::size_t n = 0; n < N; ++n) {
                C.data[m * N + n] += Aik * B.data[k * N + n];
            }
        }
    }
}


}
}
}
}
#endif // DENSEMAT_SERIAL_GEMM_IMPL_HPP
