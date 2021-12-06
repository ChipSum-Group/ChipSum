#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMV_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMV_IMPL_HPP__

#include <KokkosBlas2_gemv.hpp>


#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &y)
{

    KokkosBlas::gemv("N", 1, A, x,
                     0, y);
}


template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &y,
     const char transA[])
{

    KokkosBlas::gemv(transA, 1, A, x,
                     0, y);
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &y,
     const AlphaT& a,
     const BetaT& b
     )
{

    KokkosBlas::gemv("N", a, A, x,
                     b, y);
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::View<ValueType **> &A,
     const Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &y,
     const AlphaT& a,
     const BetaT& b,
     const char transA[])
{

    KokkosBlas::gemv(transA, a, A, x,
                     b, y);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP
