#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_AXPBY_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_AXPBY_IMPL_HPP__


#include <KokkosBlas1_axpby.hpp>


#include "../../../chipsum_macro.h"

namespace  ChipSum{

namespace  Numeric{


namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void add(
        const Kokkos::View<ValueType*>& X,
        const Kokkos::View<ValueType*>& Y) {
    KokkosBlas::axpy(1, X, Y);
}


template <typename ValueType,typename AlphaT>
CHIPSUM_FUNCTION_INLINE void axpby(
        const Kokkos::View<ValueType*>& X,
        const Kokkos::View<ValueType*>& Y,
        const AlphaT& A) {
    KokkosBlas::axpy(A, X, Y);
}

template <typename ValueType, typename AlphaT, typename BetaT>
CHIPSUM_FUNCTION_INLINE void
axpby(const Kokkos::View<ValueType*>& X,
      const Kokkos::View<ValueType*>& Y,
      const AlphaT& A,
      const BetaT& B
      ) {
    KokkosBlas::axpby(A, X, B, Y);
}
}
}
}
}


#endif // VECTOR_KOKKOSKERNELS_AXPBY_IMPL_HPP
