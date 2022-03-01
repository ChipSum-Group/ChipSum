///
/// \file     densemat_kokkoskernels_relu_impl.hpp
/// \author   Xiang Yukai
/// \group    CDCS-HPC
/// \date     
/// \brief    %stuff%
///

#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_RELU_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_RELU_IMPL_HPP__

#include <Kokkos_DualView.hpp>

#include "../chipsum/chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace  Impl {
namespace  DenseMat {

template <typename ValueType>
struct Relu_functor{
    
    Kokkos::View<ValueType **> _data;
    ValueType _alpha;

    Relu_functor(Kokkos::DualView<ValueType **> &data, ValueType alpha) : _data(data.d_view), _alpha(alpha) {}

    CHIPSUM_SPECIAL_INLINE void operator() (const int64_t i, const int64_t j) const{
        ValueType tmp = _data(i,j);
        _data(i,j) = tmp > 0 ? tmp : _alpha*tmp;
    }

};

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
relu(Kokkos::DualView<ValueType **> &A) {
     const std::size_t M = A.extent(0);
     const std::size_t N = A.extent(1);

     using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
     Relu_functor<ValueType> functor(A, 0.0);
     Kokkos::parallel_for( "relu", mdrange_policy({0,0}, {M,N}), functor);

}



template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
leakyrelu(Kokkos::DualView<ValueType **> &A) {
     const std::size_t M = A.extent(0);
     const std::size_t N = A.extent(1);

     using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
     Relu_functor<ValueType> functor(A, 0.01);
     Kokkos::parallel_for( "leakyrelu", mdrange_policy({0,0}, {M,N}), functor);

}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_RELU_IMPL_HPP
