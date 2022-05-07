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
dense(ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &input,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &weight_data,
     ChipSum::Numeric::Vector<ValueType, ChipSum::Backend::DefaultBackend> &bias_data,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &out_data) {
     
     
     auto A = input.GetData();
     auto weight = weight_data.GetData();
     auto bias = bias_data.GetData();
     auto out = out_data.GetData();

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
dense(ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &input,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &weight_data,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &out_data) {
     
     auto A = input.GetData();
     auto weight = weight_data.GetData();
     auto out = out_data.GetData();
     
     KokkosBlas::gemm("N", "N", 1, A.d_view, weight.d_view,
                     0, out.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
dense(const char transA[],
     const char transB[],
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &input,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &weight_data,
     ChipSum::Numeric::Vector<ValueType, ChipSum::Backend::DefaultBackend> &bias_data,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &out_data
     ) {

     auto A = input.GetData();
     auto weight = weight_data.GetData();
     auto bias = bias_data.GetData();
     auto out = out_data.GetData();

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
dense(const char transA[],
     const char transB[],
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &input,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &weight_data,
     ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &out_data
     ) {

     auto A = input.GetData();
     auto weight = weight_data.GetData();
     auto out = out_data.GetData();
     
     KokkosBlas::gemm(transA, transB, 1, A.d_view, weight.d_view,
                     0, out.d_view);
}

#endif // DENSEMAT_KOKKOSKERNELS_DENSE_IMPL_HPP
