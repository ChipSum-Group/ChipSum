#ifndef __CHIPSUM_TENSOR_KOKKOSKERNELS_LU_IMPL_HPP__
#define __CHIPSUM_TENSOR_KOKKOSKERNELS_LU_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>
#include "../../../chipsum_macro.h"


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Tensor {
typedef Kokkos::TeamPolicy<>               team_policy;
typedef Kokkos::TeamPolicy<>::member_type  member_type;

//Functor for Kokkos::parallel_for
template <typename ValueType>
struct LuFunctor {
    LuFunctor(const Kokkos::View<ValueType ***, Kokkos::LayoutRight> &A, const ValueType& tiny):
    _A(A),_tiny(tiny){}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const member_type &teamMember) const {
        const int i = teamMember.league_rank();
        auto A_sub = Kokkos::subview(_A, i, Kokkos::ALL(),Kokkos::ALL());
        KokkosBatched::TeamLU<member_type,KokkosBatched::Algo::LU::Unblocked>::invoke(teamMember,A_sub,_tiny);
    }

    Kokkos::View<ValueType ***, Kokkos::LayoutRight> _A;
    ValueType _tiny;
};

//Tensor3
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_lu(const Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A, const ValueType& tiny) {
    const int N = A.extent(0);
    team_policy policy(N, Kokkos::AUTO());
    LuFunctor<ValueType> functor(A.d_view,tiny);
    Kokkos::parallel_for( "LU", policy,functor);
    Kokkos::fence();
}

//Tensor4
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_lu(const Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &A, const ValueType& tiny) {
    for(size_t i=0; i<A.extent(0); ++i){
        auto A_sub = Kokkos::subview(A, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
        lu(A_sub,tiny);
    }
}

} //Tensor
}
}
}
#endif // TENSOR_KOKKOSKERNELS_LU_IMPL_HPP
