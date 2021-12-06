#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_SCAL_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_SCAL_IMPL_HPP__

#include <KokkosBlas1_scal.hpp>


#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType, typename... Args>

CHIPSUM_FUNCTION_INLINE void scal(
        const Kokkos::View<ValueType **> &A,
        Kokkos::View<ValueType **> &R,
        const ValueType alpha) {
    KokkosBlas::scal(R, alpha, A);
}

}
}
}
}

#endif // DENSEMAT_KOKKOSKERNELS_SCAL_IMPL_HPP
