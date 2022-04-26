#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_SUM_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_SUM_IMPL_HPP__

#include<typeinfo>
#include<KokkosBlas1_sum.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"
namespace ChipSum{
namespace Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE ValueType
sum(const Kokkos::DualView<ValueType *>& X) {
    return KokkosBlas::sum(X.d_view);
}


}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_SUM_IMPL_HPP
