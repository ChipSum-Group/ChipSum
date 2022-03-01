#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_DOT_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_DOT_IMPL_HPP__



#include <KokkosBlas1_dot.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE ValueType dot(
        const Kokkos::DualView<ValueType *>& x,
        const Kokkos::DualView<ValueType *>& y
        )
{

    return KokkosBlas::dot(x.d_view, y.d_view);
}

template <typename ValueType,typename RefType>
CHIPSUM_FUNCTION_INLINE void
dot(
        const Kokkos::DualView<ValueType *>& x,
        const Kokkos::DualView<ValueType *>& y,
        RefType& r
        ) {

    KokkosBlas::dot(r, x.d_view, y.d_view);
}

}
}
}
}
#endif
