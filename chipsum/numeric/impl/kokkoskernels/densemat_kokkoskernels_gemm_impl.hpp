#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP__

#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_DualView.hpp>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C) {

    KokkosBlas::gemm("N", "N", 1, A.d_view, B.d_view,
                     0, C.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C,
     const char transA[],
     const char transB[]
     ) {

    KokkosBlas::gemm(transA, transB, 1, A.d_view, B.d_view,
                     0, C.d_view);
}

template <typename ValueType,typename AlphaT,typename BetaT>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C,
     const AlphaT& a,
     const BetaT& b
     ) {

    KokkosBlas::gemm("N", "N", a, A.d_view, B.d_view,
                     b, C.d_view);
}

template <typename ValueType,typename AlphaT,typename BetaT>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C,
     const AlphaT& a,
     const BetaT& b,
     const char transA[],
     const char transB[]
     ) {

    KokkosBlas::gemm(transA, transB, a, A.d_view, B.d_view,
                     b, C.d_view);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP
