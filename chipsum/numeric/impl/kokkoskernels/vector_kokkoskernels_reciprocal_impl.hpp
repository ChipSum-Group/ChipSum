#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_RECIPROCAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_RECIPROCAL_IMPL_HPP__

#include<KokkosBlas1_reciprocal.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"
namespace ChipSum{
namespace Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
reciprocal(const Kokkos::DualView<ValueType *>& X,
           const Kokkos::DualView<ValueType *>& Y) {
    KokkosBlas::reciprocal(Y.d_view,X.d_view);
}

}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_RECIPROCAL_IMPL_HPP
