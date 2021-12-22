#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_SCAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_SCAL_IMPL_HPP__


#include <KokkosBlas1_scal.hpp>


#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace  Impl{
namespace  Vector{

template <typename ValueType>
/// y = a*x 的函数子
struct scal_functor {
    scal_functor(const Kokkos::View<ValueType*>& r,
                 const Kokkos::View<ValueType>& a,
                 const Kokkos::View<ValueType*>& x
                 )
        :__r(r),__a(a),__x(x)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
        __r(i) = __a() * __x(i);
    }

private:
    Kokkos::View<ValueType*> __r;
    Kokkos::View<ValueType> __a;
    Kokkos::View<ValueType*> __x;

};

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void scal(

        const Kokkos::View<ValueType*>& X,
        const Kokkos::View<ValueType>& A,
        const Kokkos::View<ValueType*>& R) {

    Kokkos::parallel_for(R.extent(0),
                         scal_functor<
                         ValueType
                         >(R,A,X));
}

template <typename ValueType,typename RefType>
CHIPSUM_FUNCTION_INLINE void scal(
        const Kokkos::View<ValueType*>& X,
        const RefType& a,
        const Kokkos::View<ValueType*>& R
        )
{
    KokkosBlas::scal(R, a, X);
}

}
}
}
}
#endif // VECTOR_KOKKOSKERNELS_SCAL_IMPL_HPP
