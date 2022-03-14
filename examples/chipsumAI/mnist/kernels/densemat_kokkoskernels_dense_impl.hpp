///
/// \file     densemat_kokkoskernels_dense_impl.hpp
/// \author   Xiang Yukai
/// \group    CDCS-HPC
/// \date     
/// \brief    %stuff%
///

#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_DENSE_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_DENSE_IMPL_HPP__

#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_DualView.hpp>

#include "../chipsum/chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {


template <typename ValueType>
struct Dense_functor{
    
    Kokkos::View<ValueType **> _data;
    Kokkos::View<ValueType *> _bias;

    CHIPSUM_SPECIAL_INLINE 
    Dense_functor(const Kokkos::DualView<ValueType **> &data, const Kokkos::DualView<ValueType *> &bias) : 
                    _data(data.d_view), _bias(bias.d_view) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const int64_t i, const int64_t j) const{
         _data(i,j) += _bias(j);
    }

};


template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
dense(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &weight,
     const Kokkos::DualView<ValueType *> &bias,
     Kokkos::DualView<ValueType **> &out) {

     KokkosBlas::gemm("N", "N", 1, A.d_view, weight.d_view,
                     0, out.d_view);
     
     using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
     const std::size_t M = out.extent(0);// out row
     const std::size_t N = out.extent(1);// out col

     Dense_functor<ValueType> functor(out, bias);
     Kokkos::parallel_for( "dense", mdrange_policy({0,0}, {M,N}), functor);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
dense(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &weight,
     Kokkos::DualView<ValueType **> &out) {

     KokkosBlas::gemm("N", "N", 1, A.d_view, weight.d_view,
                     0, out.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
dense(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &weight,
     const Kokkos::DualView<ValueType *> &bias,
     Kokkos::DualView<ValueType **> &out,
     const char transA[],
     const char transB[]
     ) {

     KokkosBlas::gemm(transA, transB, 1, A.d_view, weight.d_view,
                     0, out.d_view);
     
     using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
     const std::size_t M = out.extent(0);// out row
     const std::size_t N = out.extent(1);// out col
     
     Dense_functor<ValueType> functor(out, bias);
     Kokkos::parallel_for( "dense_tParameter", mdrange_policy({0,0}, {M,N}), functor);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
dense(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &weight,
     Kokkos::DualView<ValueType **> &out,
     const char transA[],
     const char transB[]
     ) {

     KokkosBlas::gemm(transA, transB, 1, A.d_view, weight.d_view,
                     0, out.d_view);
}


}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_DENSE_IMPL_HPP
