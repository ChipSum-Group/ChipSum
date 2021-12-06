#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP__

#include <KokkosBlas3_gemm.hpp>


#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType **> &B,
     Kokkos::View<ValueType **> &C) {

    KokkosBlas::gemm("N", "N", 1, A, B,
                     0, C);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType **> &B,
     Kokkos::View<ValueType **> &C,
     const char transA[],
     const char transB[]
     ) {

    KokkosBlas::gemm(transA, transB, 1, A, B,
                     0, C);
}

template <typename ValueType,typename AlphaT,typename BetaT>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType **> &B,
     Kokkos::View<ValueType **> &C,
     const AlphaT& a,
     const BetaT& b
     ) {

    KokkosBlas::gemm("N", "N", a, A, B,
                     b, C);
}

template <typename ValueType,typename AlphaT,typename BetaT>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType **> &B,
     Kokkos::View<ValueType **> &C,
     const AlphaT& a,
     const BetaT& b,
     const char transA[],
     const char transB[]
     ) {

    KokkosBlas::gemm(transA, transB, a, A, B,
                     b, C);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP
