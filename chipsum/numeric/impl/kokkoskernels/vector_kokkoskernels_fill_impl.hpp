#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_FILL_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_FILL_IMPL_HPP__

#include <KokkosBlas1_fill.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"
namespace ChipSum{
namespace Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
fill(const Kokkos::DualView<ValueType *>& X,const ValueType& alpha) {
    KokkosBlas::fill(X.d_view,alpha);
}

}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_FILL_IMPL_HPP
