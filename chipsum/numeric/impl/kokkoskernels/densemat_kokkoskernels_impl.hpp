///
/// \file     densemat_kokkoskernels_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-05
/// \brief    %stuff%
///

#ifndef __CHIPSUM_DENSEMAT_KOKKOEKERNELS_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOEKERNELS_IMPL_HPP__



#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"

#include "densemat_kokkoskernels_scal_impl.hpp"
#include "densemat_kokkoskernels_gemv_impl.hpp"
#include "densemat_kokkoskernels_gemm_impl.hpp"


/// kokkos后端由于设计失误，导致不太适合用于数值系统求解
/// 目前todo的工作是将数据结构修改为multi vector
/// 这样做的目的，举个例子，是能够例如将指定的i行抽取出来，
/// 与指定的j行进行axpby，然后方便地实现高斯消元操作

/// 如果继续按照目前的view<type**>结构走下去，
/// 在做chipsum.benchmark时势必会受到很大的影响。



static int matrix_name = 0;

namespace ChipSum {

namespace Numeric {

template <typename ValueType,typename ...Props>
struct DenseMatrix_Traits<ValueType, ChipSum::Backend::KokkosKernels,Props...>
        : public Operator_Traits<ValueType> {

    using matrix_type = Kokkos::View<ValueType **>;

    using size_type = typename matrix_type::size_type;

    using value_type = typename matrix_type::value_type;
};

namespace Impl {
namespace DenseMat {

template<typename ValueType> using traits = DenseMatrix_Traits<ValueType,ChipSum::Backend::KokkosKernels>;

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void create(
        Kokkos::View<ValueType **> &A,
        const ::std::size_t M,
        const ::std::size_t N)
{

    A = Kokkos::View<ValueType **>("densemat_" + std::to_string(matrix_name++), M,
                                N);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void create(
        Kokkos::View<ValueType **> &A,
        const std::size_t M,
        const std::size_t N,
        ValueType *src
        )
{

    if(A.extent(0)!=M||A.extent(1)!=N)
    {
        A = Kokkos::View<ValueType **>("densemat_" + ::std::to_string(matrix_name++), M,
                                    N);
    }

    typename Kokkos::View<ValueType**>::HostMirror h_A = Kokkos::create_mirror_view(A);

    using mdrange_policy = Kokkos::MDRangePolicy< Kokkos::Rank<2>,Kokkos::OpenMP > ;
      Kokkos::parallel_for( "init_A", mdrange_policy({0,0}, {M,N}), KOKKOS_LAMBDA ( const int i , const int j ) {
          h_A(i,j) = src[N*i+j];
        }
      );

    Kokkos::deep_copy(A, h_A);
}


template <typename ValueType>

CHIPSUM_FUNCTION_INLINE ValueType &
get_item(Kokkos::View<ValueType **> &A,
         const std::size_t i,
         const std::size_t j) {
    return A(i, j);
}




template <typename ValueType,typename IDT>
CHIPSUM_FUNCTION_INLINE void
set_row(const Kokkos::View<ValueType**>& A,
          const Kokkos::View<ValueType*>& a,
          const IDT&  i
          ){

    auto A_sub = Kokkos::subview(A,i,Kokkos::ALL());

    Kokkos::deep_copy(A_sub,a);

}


template <typename ValueType,typename IDT>
CHIPSUM_FUNCTION_INLINE void
set_items(const Kokkos::View<ValueType**>& A,
          const Kokkos::View<ValueType*>& a,
          const IDT&  i
          ){

    auto A_sub = Kokkos::subview(A,i,Kokkos::ALL());

    Kokkos::deep_copy(A_sub,a);

}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void print(Kokkos::View<ValueType **> &A,
                                   ::std::ostream &out) {
    ::std::size_t M = A.extent(0);
    ::std::size_t N = A.extent(1);
    auto h_A = Kokkos::create_mirror_view(A);

    Kokkos::deep_copy(h_A, A);

    cout << A.label()  <<"("
         << A.extent(0) <<","
         << A.extent(1) <<")"
         << ":" << endl;

    for (std::size_t i = 0; i < M; ++i) {
        out << " "
            << "[";
        for (std::size_t j = 0; j < N - 1; ++j) {
            out << h_A(i, j) << ", ";
        }
        out << h_A(i, N - 1) << "]" << endl;
    }
    out << endl;
}

} // namespace DenseMat
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
