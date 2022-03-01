#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMV_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMV_IMPL_HPP__

#include <KokkosBlas2_gemv.hpp>
#include <Kokkos_DualView.hpp>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y)
{

    KokkosBlas::gemv("N", 1, A.d_view, x.d_view,
                     0, y.d_view);
}


template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y,
     const char transA[])
{

    KokkosBlas::gemv(transA, 1, A.d_view, x.d_view,
                     0, y.d_view);
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y,
     const AlphaT& a,
     const BetaT& b
     )
{

    KokkosBlas::gemv("N", a, A.d_view, x.d_view,
                     b, y.d_view);
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y,
     const AlphaT& a,
     const BetaT& b,
     const char transA[])
{

    KokkosBlas::gemv(transA, a, A.d_view, x.d_view,
                     b, y.d_view);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP
