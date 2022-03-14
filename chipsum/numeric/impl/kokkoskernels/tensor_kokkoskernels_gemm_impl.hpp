#ifndef __CHIPSUM_TENSOR_KOKKOSKERNELS_GEMM_IMPL_HPP__
#define __CHIPSUM_TENSOR_KOKKOSKERNELS_GEMM_IMPL_HPP__

#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_DualView.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Handle.hpp>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Tensor {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_gemm(const Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A,
     const Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &B,
     Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &C) {

     KokkosBatched::BatchedGemmHandle batchedGemmHandle;
     
     KokkosBatched::BatchedGemm<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose, KokkosBatched::BatchLayout::Left>
                              (&batchedGemmHandle, 1, A.d_view, B.d_view, 0, C.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_gemm(const Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &A,
     const Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &B,
     Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &C) {

     KokkosBatched::BatchedGemmHandle batchedGemmHandle;
     
     for(size_t i=0; i<A.extent(0); ++i){
          
          auto A_sub = Kokkos::subview(A.d_view, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          auto B_sub = Kokkos::subview(B.d_view, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          auto C_sub = Kokkos::subview(C.d_view, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          KokkosBatched::BatchedGemm<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose, KokkosBatched::BatchLayout::Left>
                              (&batchedGemmHandle, 1, A_sub, B_sub, 0, C_sub);
     }
}

}
}
}
}
#endif // TENSOR_KOKKOSKERNELS_GEMM_IMPL_HPP
