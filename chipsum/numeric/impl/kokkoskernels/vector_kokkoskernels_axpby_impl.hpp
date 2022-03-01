#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_AXPBY_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_AXPBY_IMPL_HPP__


#include <KokkosBlas1_axpby.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"

namespace  ChipSum{

namespace  Numeric{


namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void add(
        const Kokkos::DualView<ValueType *>& X,
        Kokkos::DualView<ValueType *>& Y) {
    KokkosBlas::axpy(1, X.d_view, Y.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void axpby(
        const Kokkos::DualView<ValueType *>& X,
        Kokkos::DualView<ValueType *>& Y) {
    KokkosBlas::axpy(1, X.d_view, Y.d_view);
}


template <typename ValueType,typename AlphaT>
CHIPSUM_FUNCTION_INLINE void axpby(
        const Kokkos::DualView<ValueType *>& X,
        Kokkos::DualView<ValueType *>& Y,
        const AlphaT& A) {
    KokkosBlas::axpy(A, X.d_view, Y.d_view);
}

template <typename ValueType, typename AlphaT, typename BetaT>
CHIPSUM_FUNCTION_INLINE void
axpby(const Kokkos::DualView<ValueType *>& X,
      Kokkos::DualView<ValueType *>& Y,
      const AlphaT& A,
      const BetaT& B
      ) {
    KokkosBlas::axpby(A, X.d_view, B, Y.d_view);
}
}
}
}
}


#endif // VECTOR_KOKKOSKERNELS_AXPBY_IMPL_HPP
