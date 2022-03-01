#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_LU_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_LU_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include <KokkosBatched_LU_Decl.hpp>
// #include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>

#include "../../../chipsum_macro.h"
// using namespace KokkosBatched;


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
lu(const Kokkos::DualView<ValueType **> &A, const ValueType& tiny) {
    KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Unblocked>::invoke(A.h_view,tiny);
    Kokkos::deep_copy(A.d_view, A.h_view);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_LU_IMPL_HPP
