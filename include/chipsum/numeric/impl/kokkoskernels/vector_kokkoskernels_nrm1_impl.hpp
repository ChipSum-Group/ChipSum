#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_NRM1_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_NRM1_IMPL_HPP__

#include <KokkosBlas1_nrm1.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE auto
norm1(const Kokkos::DualView<ValueType *>& X) {
    return KokkosBlas::nrm1(X.d_view);
}

template <typename ValueType,typename RefType>
CHIPSUM_FUNCTION_INLINE auto
norm1(const Kokkos::DualView<ValueType *>& X,RefType& acc) {
    return KokkosBlas::nrm1(acc,X.d_view);
}
}
}
}
}

#endif // VECTOR_KOKKOSKERNELS_NRM1_IMPL_HPP
