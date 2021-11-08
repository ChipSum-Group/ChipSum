///
/// \file     densemat_kokkoskernels_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-05
/// \brief    %stuff%
///

#ifndef __CHIPSUM_DENSEMAT_KOKKOEKERNELS_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOEKERNELS_IMPL_HPP__

#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas3_gemm.hpp>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"


/// kokkos后端由于设计失误，导致不太适合用于数值系统求解
/// 目前todo的工作是将数据结构修改为multi vector
/// 这样做的目的，举个例子，是能够例如将指定的i行抽取出来，
/// 与指定的j行进行axpby，然后方便地实现高斯消元操作

/// 如果继续按照目前的view<type**>结构走下去，
/// 在做chipsum.benchmark时势必会受到很大的影响。



static int matrix_name = 0;

namespace ChipSum {

namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct DenseMatrix_Traits<ScalarType, SizeType, ChipSum::Backend::KokkosKernels,
                          Props...>
    : public Operator_Traits<ScalarType, SizeType,
                             ChipSum::Backend::KokkosKernels, Props...> {

  using matrix_type = Kokkos::View<double **>;

  using size_type = std::size_t;
};

namespace Impl {
namespace DenseMat {

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void create(const std::size_t M, const std::size_t N,
                                    Kokkos::View<double **> &A) {

  A = Kokkos::View<double **>("densemat_" + std::to_string(matrix_name++), M,
                                N);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void fill(const std::size_t M, const std::size_t N,
                                  ScalarType *src,
                                  Kokkos::View<double **> &dst) {
  auto h_dst = Kokkos::View<double **>(src, M, N);

  Kokkos::deep_copy(dst, h_dst);
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void
mult(const std::size_t M, const std::size_t N, const std::size_t K,
     const Kokkos::View<double **> &A, const Kokkos::View<double **> &B,
     Kokkos::View<double **> &C) {

  KokkosBlas::gemm("N", "N", static_cast<ScalarType>(1), A, B,
                   static_cast<ScalarType>(0), C);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
mult(const std::size_t M, const std::size_t N, const Kokkos::View<double **> &A,
     const Kokkos::View<double *> &x, Kokkos::View<double *> &y) {

  KokkosBlas::gemv("N", static_cast<ScalarType>(1), A, x,
                   static_cast<ScalarType>(0), y);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void scal(const ScalarType alpha, const std::size_t M,
                                  const std::size_t N,
                                  Kokkos::View<double **> &A) {
  KokkosBlas::scal(A, alpha, A);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType &
get_item(const std::size_t i, const std::size_t j, const std::size_t M,
        const std::size_t N, Kokkos::View<double **> &A) {
  return A(i, j);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void print(const std::size_t M, const std::size_t N,
                                   Kokkos::View<double **> &A,
                                   std::ostream &out) {
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
    out << h_A(i, M - 1) << "]" << endl;
  }
  out << endl;
}

} // namespace DenseMat
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
