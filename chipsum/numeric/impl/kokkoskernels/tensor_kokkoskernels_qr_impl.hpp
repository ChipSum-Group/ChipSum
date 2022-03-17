#ifndef __CHIPSUM_TENSOR_KOKKOSKERNELS_QR_IMPL_HPP__
#define __CHIPSUM_TENSOR_KOKKOSKERNELS_QR_IMPL_HPP__

#include <Kokkos_DualView.hpp>
#include <KokkosBatched_QR_Decl.hpp>
#include <KokkosBatched_QR_Serial_Impl.hpp>
#include <KokkosBatched_QR_TeamVector_Impl.hpp>
#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Tensor {

typedef Kokkos::TeamPolicy<>               team_policy;
typedef Kokkos::TeamPolicy<>::member_type  member_type;

//Functor for Kokkos::parallel_for
template <typename ValueType>
struct QrFunctor {
    QrFunctor(const Kokkos::View<ValueType ***, Kokkos::LayoutRight> &A, 
              const Kokkos::View<ValueType **>  &x,
              const Kokkos::View<ValueType **>  &y)
              :_A(A),_x(x),_y(y){}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const member_type &teamMember) const {
        const int i = teamMember.league_rank();
        auto A_sub = Kokkos::subview(_A, i, Kokkos::ALL(),Kokkos::ALL());
        auto x_sub = Kokkos::subview(_x, i, Kokkos::ALL());
        auto y_sub = Kokkos::subview(_y, i, Kokkos::ALL());
        KokkosBatched::TeamVectorQR<member_type,KokkosBatched::Algo::QR::Unblocked>
        ::invoke(teamMember,A_sub,x_sub,y_sub);
    }

    Kokkos::View<ValueType ***, Kokkos::LayoutRight> _A;
    Kokkos::View<ValueType **> _x;
    Kokkos::View<ValueType **> _y;
};

//Tensor3
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_qr( Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A, 
          Kokkos::DualView<ValueType **>  &x,
          Kokkos::DualView<ValueType **>  &y) {
    const int N = A.extent(0);
    team_policy policy(N,Kokkos::AUTO());
    QrFunctor<ValueType> functor(A.d_view,x.d_view,y.d_view);
    Kokkos::parallel_for( "QR", policy,functor);
    Kokkos::fence();
}

//Tensor4
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
batch_qr(Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &A,
         Kokkos::DualView<ValueType ***>  &x,
         Kokkos::DualView<ValueType ***>  &y) {
    for(size_t i=0; i<A.extent(0); ++i){
        auto A_sub = Kokkos::subview(A, i, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
        auto x_sub = Kokkos::subview(x, i, Kokkos::ALL(),Kokkos::ALL());
        auto y_sub = Kokkos::subview(y, i, Kokkos::ALL(),Kokkos::ALL());
        batch_qr(A_sub,x,y);
    }
}


}//Tensor
}
}
}
#endif //TENSOR_KOKKOSKERNELS_QR_IMPL_HPP
