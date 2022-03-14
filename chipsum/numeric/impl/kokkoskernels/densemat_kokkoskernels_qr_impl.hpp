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
// SERIAL CODE
// template <typename ValueType>
// CHIPSUM_FUNCTION_INLINE void
// qr(const Kokkos::DualView<ValueType **> &A, 
//    const Kokkos::DualView<ValueType *> &x,
//    const Kokkos::DualView<ValueType *> &y) {
//     KokkosBatched::SerialQR<KokkosBatched::Algo::QR::Unblocked>::invoke(A.h_view,x.h_view,y.h_view);
//     Kokkos::deep_copy(A.d_view, A.h_view);
//     Kokkos::deep_copy(x.d_view, x.h_view);
//     Kokkos::deep_copy(y.d_view, y.h_view);
// }

typedef Kokkos::TeamPolicy<>               team_policy;
typedef Kokkos::TeamPolicy<>::member_type  member_type;

template <typename ValueType>
struct Parallel_QR {
    Parallel_QR(const Kokkos::DualView<ValueType **> &A, 
                const Kokkos::DualView<ValueType *> &x,
                const Kokkos::DualView<ValueType *> &y):_A(A),_x(x),_y(y){}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const member_type &teamMember) const {
        KokkosBatched::TeamVectorQR<member_type,KokkosBatched::Algo::QR::Unblocked>::invoke(teamMember,_A.d_view,_x.d_view,_y.d_view);
    }

    Kokkos::DualView<ValueType **> _A;
    Kokkos::DualView<ValueType *> _x;
    Kokkos::DualView<ValueType *> _y;
};

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
qr(const Kokkos::DualView<ValueType **> &A, 
   const Kokkos::DualView<ValueType *> &x,
   const Kokkos::DualView<ValueType *> &y) {
    int N = 1;
    team_policy policy(N,Kokkos::AUTO());
    Parallel_QR<ValueType> functor(A,x,y);
    Kokkos::parallel_for( "QR", policy,functor);
    Kokkos::fence();
}







}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_QR_IMPL_HPP
