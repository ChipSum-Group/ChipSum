#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IAMAX_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IAMAX_IMPL_HPP__

#include<typeinfo>
#include<KokkosBlas1_iamax.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"
namespace ChipSum{
namespace Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE
typename Kokkos::View<ValueType *>::size_type
iamax(const Kokkos::DualView<ValueType *>& X) {
    return KokkosBlas::iamax(X.d_view);
}

template <typename ValueType,typename SizeType>
CHIPSUM_FUNCTION_INLINE void
iamax(const Kokkos::DualView<ValueType *>& X,
      const typename Kokkos::View<SizeType>& r) {
    KokkosBlas::iamax(r, X.d_view);
}



}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_IAMAX_IMPL_HPP
