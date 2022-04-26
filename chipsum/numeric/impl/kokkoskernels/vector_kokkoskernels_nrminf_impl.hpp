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

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
norminf(const Kokkos::DualView<ValueType *>& X,
        const typename Kokkos::View<ValueType>& r) {
    // KokkosBlas::nrminf(r, X.d_view);
    std::cout<<"todo ... "<<std::endl;
}



}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_NRMINF_IMPL_HPP
