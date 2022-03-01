///
/// \file     densemat_kokkoskernels_norm_impl.hpp
/// \author   Xiang Yukai
/// \group    CDCS-HPC
/// \date     
/// \brief    %stuff%
///

#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_NORM_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_NORM_IMPL_HPP__

#include <Kokkos_DualView.hpp>

#include "../chipsum/chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace  Impl {
namespace  DenseMat {

template <typename ValueType>
struct Normsum_functor{
    
    Kokkos::View<ValueType **> _data;

    CHIPSUM_SPECIAL_INLINE 
    Normsum_functor(Kokkos::DualView<ValueType **> &data) : _data(data.d_view) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const int64_t i, const int64_t j, ValueType & tmpSum) const{
        tmpSum += _data(i,j);
    }

};

template <typename ValueType>
struct Norm_functor{
    
    Kokkos::View<ValueType **> _data;
    ValueType _sum;

    CHIPSUM_SPECIAL_INLINE 
    Norm_functor(Kokkos::DualView<ValueType **> &data, ValueType &sum) : _data(data.d_view), _sum(sum) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const int64_t i, const int64_t j) const{
        // printf("%f*_sum* ", _sum);// check reduce和是否正确
        _data(i,j) /= _sum;
    }

};

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
norm(Kokkos::DualView<ValueType **>& A) {
    const std::size_t M = A.extent(0);// row
    const std::size_t N = A.extent(1);// col

    ValueType val_sum = 0.0;
    using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Normsum_functor<ValueType> functor_sum(A);
    Kokkos::parallel_reduce( "norm_sum", mdrange_policy({0,0}, {M,N}), functor_sum, val_sum);

    Norm_functor<ValueType> functor(A, val_sum);
    Kokkos::parallel_for( "norm", mdrange_policy({0,0}, {M,N}), functor);
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_NORM_IMPL_HPP
