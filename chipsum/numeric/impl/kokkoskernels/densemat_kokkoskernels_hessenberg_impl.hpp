#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_HESSENBERG_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_HESSENBERG_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include <KokkosBatched_Hessenberg_Serial_Internal.hpp>
#include "../../../chipsum_macro.h"
// using namespace KokkosBatched;


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
hessenberg(const Kokkos::DualView<ValueType **> &A, 
    const Kokkos::DualView<ValueType *> &t,
    const Kokkos::DualView<ValueType *> &w) {
    KokkosBatched::SerialHessenbergInternal::invoke(A.h_view.extent(0), A.h_view.extent(1), 
        A.h_view.data(), A.h_view.stride_0(), A.h_view.stride_1(),
        t.h_view.data(), t.h_view.stride_0(), w.h_view.data());
    // std::cout<<"A.stride_0: "<<A.h_view.stride_0()<<" A.stride_1: "<<A.h_view.stride_1()<<std::endl;
    // std::cout<<"t.stride_0: "<<t.h_view.stride_0()<<std::endl;
    Kokkos::deep_copy(A.d_view, A.h_view);
    Kokkos::deep_copy(t.d_view, t.h_view);
    Kokkos::deep_copy(w.d_view, w.h_view);

}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_HESSENBERG_IMPL_HPP
