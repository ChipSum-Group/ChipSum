#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_TRMM_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_TRMM_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include<KokkosBlas3_trmm.hpp>
#include "../../../chipsum_macro.h"


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
trmm(Kokkos::DualView<ValueType **> &A,
     Kokkos::DualView<ValueType **> &B, 
     const ValueType &alpha,
     const char side[],
     const char uplo[],
     const char trans[],
     const char diag[]) {
    KokkosBlas::trmm(side,uplo,trans,diag,alpha,A.d_view,B.d_view);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_TRMM_IMPL_HPP
