#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_LU_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_LU_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>

#include "../../../chipsum_macro.h"
// using namespace KokkosBatched;


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {
//serial code
// template <typename ValueType>
// CHIPSUM_FUNCTION_INLINE void
// lu(const Kokkos::DualView<ValueType **> &A, const ValueType& tiny) {
//     KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Unblocked>::invoke(A.h_view,tiny);
//     Kokkos::deep_copy(A.d_view, A.h_view);
// }

typedef Kokkos::TeamPolicy<>               team_policy;
typedef Kokkos::TeamPolicy<>::member_type  member_type;

template <typename ValueType>
struct LuFunctor {
    LuFunctor(const Kokkos::View<ValueType **> &A, const ValueType& tiny):_A(A),_tiny(tiny){}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const member_type &teamMember) const {
        KokkosBatched::TeamLU<member_type,KokkosBatched::Algo::LU::Unblocked>::invoke(teamMember,_A,_tiny);
    }

    Kokkos::View<ValueType **> _A;
    ValueType _tiny;
};

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
lu(const Kokkos::DualView<ValueType **> &A, const ValueType& tiny) {
    const int N = 1;
    team_policy policy(N, Kokkos::AUTO());
    LuFunctor<ValueType> functor(A.d_view,tiny);
    Kokkos::parallel_for( "LU", policy,functor);
    Kokkos::fence();
}


}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_LU_IMPL_HPP
