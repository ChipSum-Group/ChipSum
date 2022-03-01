#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_QR_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_QR_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include <KokkosBatched_QR_Decl.hpp>
#include <KokkosBatched_QR_Serial_Impl.hpp>
#include <KokkosBatched_QR_TeamVector_Impl.hpp>
#include "../../../chipsum_macro.h"
// using namespace KokkosBatched;


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
qr(const Kokkos::DualView<ValueType **> &A, 
   const Kokkos::DualView<ValueType *> &x,
   const Kokkos::DualView<ValueType *> &y) {
    KokkosBatched::SerialQR<KokkosBatched::Algo::QR::Unblocked>::invoke(A.h_view,x.h_view,y.h_view);
    Kokkos::deep_copy(A.d_view, A.h_view);
    Kokkos::deep_copy(x.d_view, x.h_view);
    Kokkos::deep_copy(y.d_view, y.h_view);

}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_QR_IMPL_HPP
