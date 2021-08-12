/*
 * @Description: 稠密矩阵dense_matrix的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:42:26
 */

#ifndef __CHIPSUM_DENSEMAT_KOKKOEKERNELS_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOEKERNELS_IMPL_HPP__

#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas3_gemm.hpp>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

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
/**
 * @brief Create
 * @param M
 * @param N
 * @param mat
 */
CHIPSUM_FUNCTION_INLINE void Create(const std::size_t M, const std::size_t N,
                                    Kokkos::View<double **> &mat) {

  mat = Kokkos::View<double **>("densemat_" + std::to_string(matrix_name++), M,
                                N);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Fill
 * @param M
 * @param N
 * @param src
 * @param dst
 */
CHIPSUM_FUNCTION_INLINE void Fill(const std::size_t M, const std::size_t N,
                                  ScalarType *src,
                                  Kokkos::View<double **> &dst) {
  auto h_dst = Kokkos::View<double **>(src, M, N);

  Kokkos::deep_copy(dst, h_dst);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Mult
 * @param M
 * @param N
 * @param K
 * @param A
 * @param B
 * @param C
 */
CHIPSUM_FUNCTION_INLINE void
Mult(const std::size_t M, const std::size_t N, const std::size_t K,
     const Kokkos::View<double **> &A, const Kokkos::View<double **> &B,
     Kokkos::View<double **> &C) {

  KokkosBlas::gemm("N", "N", 1.0, A, B, 0.0, C);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Mult
 * @param M
 * @param N
 * @param A
 * @param x
 * @param y
 */
CHIPSUM_FUNCTION_INLINE void
Mult(const std::size_t M, const std::size_t N, const Kokkos::View<double **> &A,
     const Kokkos::View<double *> &x, Kokkos::View<double *> &y) {

  KokkosBlas::gemv("N", 1.0, A, x, 0.0, y);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Scal
 * @param alpha
 * @param M
 * @param N
 * @param mat
 */
CHIPSUM_FUNCTION_INLINE void Scal(ScalarType alpha, const std::size_t M,
                                  const std::size_t N,
                                  Kokkos::View<double **> &mat) {
  KokkosBlas::scal(mat, alpha, mat);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief GetItem
 * @param i
 * @param j
 * @param M
 * @param N
 * @param mat
 * @return
 */
CHIPSUM_FUNCTION_INLINE ScalarType &
GetItem(const std::size_t i, const std::size_t j, const std::size_t M,
        const std::size_t N, Kokkos::View<double **> &mat) {
  return mat(i, j);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Print
 * @param M
 * @param N
 * @param mat
 * @param out
 */
CHIPSUM_FUNCTION_INLINE void Print(const std::size_t M, const std::size_t N,
                                   Kokkos::View<double **> &mat,
                                   std::ostream &out) {
  auto h_mat = Kokkos::create_mirror_view(mat);

  Kokkos::deep_copy(h_mat, mat);

  cout << mat.label() << ":" << endl;

  for (std::size_t i = 0; i < M; ++i) {
    out << " "
        << "[";
    for (std::size_t j = 0; j < M - 1; ++j) {
      out << h_mat(i, j) << ", ";
    }
    out << h_mat(i, M - 1) << "]" << endl;
  }
  out << endl;
}

} // namespace DenseMat
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
