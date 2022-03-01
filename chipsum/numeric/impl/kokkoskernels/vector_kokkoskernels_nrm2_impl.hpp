#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_NRM2_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_NRM2_IMPL_HPP__

#include <KokkosBlas1_nrm2.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE auto
norm2(const Kokkos::DualView<ValueType *>& X) {
    return KokkosBlas::nrm2(X.d_view);
}

}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_NRM2_IMPL_HPP
