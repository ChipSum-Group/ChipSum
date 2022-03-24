#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_NRMINF_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_NRMINF_IMPL_HPP__

#include <KokkosBlas1_nrminf.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"
namespace ChipSum{
namespace Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE auto
norminf(const Kokkos::DualView<ValueType *>& X) {
    return KokkosBlas::nrminf(X.d_view);
}

}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_NRMINF_IMPL_HPP
