#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_TRTRI_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_TRTRI_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include<KokkosBlas_trtri.hpp>
#include "../../../chipsum_macro.h"


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE int
trtri(Kokkos::DualView<ValueType **> &A,
     const char uplo[],
     const char diag[]) {
    return KokkosBlas::trtri(uplo,diag,A.d_view);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_TRTRI_IMPL_HPP
