///
/// \file     tensor_kokkoskernels_impl.hpp
/// \author   Yukai Xiang
/// \group    CDCS-HPC
/// \date     2022-03-07
/// \brief    %stuff%
///

#ifndef __CHIPSUM_TENSOR_KOKKOEKERNELS_IMPL_HPP__
#define __CHIPSUM_TENSOR_KOKKOEKERNELS_IMPL_HPP__

#include <Kokkos_DualView.hpp>

#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"

// #include "densemat_kokkoskernels_scal_impl.hpp"
// #include "densemat_kokkoskernels_gemv_impl.hpp"
#include "tensor_kokkoskernels_gemm_impl.hpp"
#include "tensor_kokkoskernels_gemv_impl.hpp"
#include "tensor_kokkoskernels_lu_impl.hpp"
#include "tensor_kokkoskernels_qr_impl.hpp"

// #include "../examples/chipsumAI/mnist/kernels/densemat_kokkoskernels_relu_impl.hpp"
// #include "../examples/chipsumAI/mnist/kernels/densemat_kokkoskernels_activation_impl.hpp"
// #include "../examples/chipsumAI/mnist/kernels/densemat_kokkoskernels_norm_impl.hpp"
// #include "../examples/chipsumAI/mnist/kernels/densemat_kokkoskernels_dense_impl.hpp"

// #include "densemat_kokkoskernels_lu_impl.hpp"
// #include "densemat_kokkoskernels_qr_impl.hpp"
// #include "densemat_kokkoskernels_hessenberg_impl.hpp"



static int tensor_name = 0;

namespace ChipSum {
namespace Numeric {

template <size_t NDIM, typename ValueType>
struct Tensor_Traits_t{
    using type = typename Tensor_Traits_t<NDIM-1, ValueType>::type*;
};


template <typename ValueType>
struct Tensor_Traits_t<0, ValueType>{
    using type = ValueType;
};


template <size_t NDIM, typename ValueType, typename ...Props>
struct Tensor_Traits<NDIM, ValueType, ChipSum::Backend::KokkosKernels, Props...>
        : public Operator_Traits<ValueType> {

    using tensor_type = typename Kokkos::DualView<typename Tensor_Traits_t<NDIM, ValueType>::type, Kokkos::LayoutRight>;
    using backend_type = ChipSum::Backend::KokkosKernels;
    using size_type = typename tensor_type::size_type;
    using value_type = typename tensor_type::value_type;

};


namespace Impl {
namespace Tensor {

    // template<typename ValueType> using traits = Tensor_Traits<ValueType>;

    template<typename T, typename ...Args>
    CHIPSUM_FUNCTION_INLINE void create(Kokkos::DualView<T, Kokkos::LayoutRight> & dv, Args ...args){
        static_assert(dv.rank==sizeof ...(args), "Tensor shape is not match size of args");
        dv = Kokkos::DualView<T, Kokkos::LayoutRight>("tensor_" + std::to_string(tensor_name++), args...);
    }


    // 可变模板参数 给定值初始化
    template <typename T, typename ValueType, typename ...Args>
    CHIPSUM_FUNCTION_INLINE void create(
        ValueType *src,
        Kokkos::DualView<T, Kokkos::LayoutRight> &A,
        Args ...args)
    { 
        static_assert(A.rank==sizeof ...(args), "Tensor shape is not match size of args");
        A = Kokkos::DualView<T, Kokkos::LayoutRight>("tensor_" + std::to_string(tensor_name++), args...);

        Kokkos::View<T, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > tmp_view (src, args...);
        Kokkos::deep_copy(A.h_view, tmp_view);
        Kokkos::deep_copy(A.d_view, A.h_view);
    }
    // 

    /* template <typename ValueType>
    CHIPSUM_FUNCTION_INLINE void create(
        ValueType *src,
        Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A,
        const ::std::size_t B,
        const ::std::size_t M,
        const ::std::size_t N)
    {
        if(A.extent(0)!=B||A.extent(1)!=M||A.extent(2)!=N)
        {
            A = Kokkos::DualView<ValueType ***, Kokkos::LayoutRight>("tensor_" + std::to_string(tensor_name++), B, M,
                                    N);
        }
        using mdrange_policy = Kokkos::MDRangePolicy< Kokkos::Rank<3>, Kokkos::OpenMP >;
        Kokkos::parallel_for( "init_A", mdrange_policy({0,0,0}, {B,M,N}), KOKKOS_LAMBDA ( const int i , const int j, const int k ) {
            A.h_view(i,j,k) = src[i*M*N+j*N+k];
        }
        );

        Kokkos::deep_copy(A.d_view, A.h_view);
    }


    template <typename ValueType>
    CHIPSUM_FUNCTION_INLINE void create(
            ValueType *src,
            Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &A,
            const ::std::size_t K,
            const ::std::size_t B,
            const ::std::size_t M,
            const ::std::size_t N)
    {
        if(A.extent(0)!=K||A.extent(1)!=B||A.extent(2)!=M||A.extent(3)!=N)
        {
            A = Kokkos::DualView<ValueType ****, Kokkos::LayoutRight>("tensor_" + std::to_string(tensor_name++), K,
                                        B, M, N);
        }

        using mdrange_policy = Kokkos::MDRangePolicy< Kokkos::Rank<4>, Kokkos::OpenMP >;
        Kokkos::parallel_for( "init_A", mdrange_policy({0,0,0,0}, {K,B,M,N}), KOKKOS_LAMBDA ( const int i , const int j, const int k, const int p ) {
            A.h_view(i,j,k,p) = src[i*B*M*N+j*M*N+k*N+p];
        }
        );

        Kokkos::deep_copy(A.d_view, A.h_view);
    } */


    template <typename T>
    CHIPSUM_FUNCTION_INLINE void
    device_to_host(Kokkos::DualView<T, Kokkos::LayoutRight> &A) 
    {
        Kokkos::deep_copy(A.h_view, A.d_view);
    }


    template <typename T>
    CHIPSUM_FUNCTION_INLINE void
    host_to_device(Kokkos::DualView<T, Kokkos::LayoutRight> &A) 
    {
        Kokkos::deep_copy(A.d_view, A.h_view);
    }


    template <typename ValueType>
    CHIPSUM_FUNCTION_INLINE void print(Kokkos::DualView<ValueType ***, Kokkos::LayoutRight> &A,
                                    ::std::ostream &out) {
        ::std::size_t B = A.extent(0);
        ::std::size_t M = A.extent(1);
        ::std::size_t N = A.extent(2);

        Kokkos::deep_copy(A.h_view, A.d_view);

        cout << A.h_view.label()  <<"("
            << A.extent(0) <<","
            << A.extent(1) <<","
            << A.extent(2) <<")"
            << ":" << endl;

        out << " "
            << "[";
        for (std::size_t i = 0; i < B ; ++i) {
            out << " "
                << "[";
            for (std::size_t j = 0; j < M ; ++j) {
                out << " "
                    << "[";
                for (std::size_t k = 0; k < N-1 ; ++k){
                    // out << " "
                    //     << "[";
                    out << A.h_view(i, j, k) << ", ";
                }
                out << A.h_view(i, j, N - 1) << "]" << endl;
                // out << "]" << endl;
            }
            // out << A.h_view(i, M-1, N - 1) << "]" << endl;
            out << "]" << endl;
        }
        out << "]" << endl;
        out << endl;
    }


    template <typename ValueType>
    CHIPSUM_FUNCTION_INLINE void print(Kokkos::DualView<ValueType ****, Kokkos::LayoutRight> &A,
                                    ::std::ostream &out) {
        ::std::size_t K = A.extent(0);
        ::std::size_t B = A.extent(1);
        ::std::size_t M = A.extent(2);
        ::std::size_t N = A.extent(3);

        Kokkos::deep_copy(A.h_view, A.d_view);

        cout << A.h_view.label()  <<"("
            << A.extent(0) <<","
            << A.extent(1) <<","
            << A.extent(2) <<","
            << A.extent(3) <<")"
            << ":" << endl;

        out << " "
            << "[";
        for (std::size_t i = 0; i < K ; ++i) {
            for (std::size_t j = 0; j < B ; ++j) {
                for (std::size_t k = 0; k < M ; ++k){
                    for (std::size_t p = 0; p < N-1 ; ++p){
                        out << A.h_view(i, j, k, p) << " ";
                    }
                    if (k<M-1)
                        out << A.h_view(i, j, k, N - 1) << ", ";
                    else
                        out << A.h_view(i, j, k, N - 1);
                }
                out << "; ";
            }
        }
        out << "]" << endl;
        out << endl;
    }

} // namespace Tensor
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_TENSOR_BLAS_IMPL_HPP__
