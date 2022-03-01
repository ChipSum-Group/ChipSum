#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_SCAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_SCAL_IMPL_HPP__


#include <KokkosBlas1_scal.hpp>
#include <Kokkos_DualView.hpp>
#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace  Impl{
namespace  Vector{

template <typename ValueType>
/// y = a*x 的函数子
struct scal_functor {
    scal_functor(const Kokkos::DualView<ValueType *>& r,
                 const Kokkos::View<ValueType>& a,
                 const Kokkos::DualView<ValueType *>& x
                 )
        :__r(r),__a(a),__x(x)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
        __r(i) = __a() * __x(i);
    }

private:
    Kokkos::DualView<ValueType *> __r;
    Kokkos::View<ValueType> __a;
    Kokkos::DualView<ValueType *> __x;

};

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void scal(

        const Kokkos::DualView<ValueType *>& X,
        const Kokkos::View<ValueType>& A,
        const Kokkos::DualView<ValueType *>& R) {

    Kokkos::parallel_for(R.extent(0),
                         scal_functor<
                         ValueType
                         >(R,A,X));
}

template <typename ValueType,typename RefType>
CHIPSUM_FUNCTION_INLINE void scal(
        const Kokkos::DualView<ValueType *>& X,
        const RefType& a,
        const Kokkos::DualView<ValueType *>& R
        )
{
    KokkosBlas::scal(R.d_view, a, X.d_view);
}

}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_SCAL_IMPL_HPP
