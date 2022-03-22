#ifndef __CHIPSUM_TENSOR_KOKKOSKERNELS_GEMV_IMPL_HPP__
#define __CHIPSUM_TENSOR_KOKKOSKERNELS_GEMV_IMPL_HPP__

#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_DualView.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Tensor {

template <typename ValueType>
struct GemvFunctor{

    typedef Kokkos::TeamPolicy<>::member_type  member_type;
    Kokkos::View<ValueType ***, Kokkos::LayoutRight> _A;
    Kokkos::View<ValueType ***, Kokkos::LayoutRight> _B;
    Kokkos::View<ValueType ***, Kokkos::LayoutRight> _C;
    ValueType _alpha;
    ValueType _beta;

    GemvFunctor(Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A, Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &B, 
               Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &C, ValueType alpha, ValueType beta) : 
               _A(A.d_view), _B(B.d_view), _C(C.d_view), _alpha(alpha), _beta(beta) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const member_type & teamMember) const{

          const int i = teamMember.league_rank();
        
          auto A_sub = Kokkos::subview(_A, i, Kokkos::ALL(),Kokkos::ALL());
          auto B_sub = Kokkos::subview(_B, i, Kokkos::ALL(),0);
          auto C_sub = Kokkos::subview(_C, i, Kokkos::ALL(),0);

          KokkosBatched::TeamGemv<member_type, KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Gemv::Blocked>::invoke(teamMember, _alpha, A_sub, B_sub, _beta, C_sub);
    }
};



template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_gemv(Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A,
     Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &B,
     Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &C) {
     
     typedef Kokkos::TeamPolicy<> team_policy;

     team_policy policy(A.d_view.extent(0), Kokkos::AUTO);
     GemvFunctor<ValueType> functor(A, B, C, 1.0, 0.0);

     Kokkos::parallel_for( "batch_gemv", policy, functor);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_gemv(Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &A,
     Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &B,
     Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &C) {
     
     for(size_t i=0; i<A.extent(0); ++i){
          
          // auto A_sub = Kokkos::subview(A, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          // auto B_sub = Kokkos::subview(B, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          // auto C_sub = Kokkos::subview(C, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

          // auto A_sub = Kokkos::subview(A, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          // auto B_sub = Kokkos::subview(B, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
          // auto C_sub = Kokkos::subview(C, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

          Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> A_sub(A, i, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
          Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> B_sub(B, i, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
          Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> C_sub(C, i, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
          
          typedef Kokkos::TeamPolicy<> team_policy;

          team_policy policy(A_sub.d_view.extent(0), Kokkos::AUTO);
          GemvFunctor<ValueType> functor(A_sub, B_sub, C_sub, 1.0, 0.0);

          Kokkos::parallel_for( "batch_gemv_4d", policy, functor);
     }
}

}
}
}
}
#endif // TENSOR_KOKKOSKERNELS_GEMV_IMPL_HPP
