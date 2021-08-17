/*
 * @Description: 稠密矩阵dense_matrix的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-16 15:59:57
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
 * @description: 创建未初始化的稠密矩阵
 * @param {*} M 行数
 * @param {*} N 列数
 * @param {*} A 稠密矩阵（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Create(const std::size_t M, const std::size_t N,
                                    Kokkos::View<double **> &A) {

  A = Kokkos::View<double **>("densemat_" + std::to_string(matrix_name++), M,
                                N);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建初始化的稠密矩阵
 * @param {*} M 行数
 * @param {*} N 列数
 * @param {*} src POD数据源
 * @param {*} dst 稠密矩阵（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Fill(const std::size_t M, const std::size_t N,
                                  ScalarType *src,
                                  Kokkos::View<double **> &dst) {
  auto h_dst = Kokkos::View<double **>(src, M, N);

  Kokkos::deep_copy(dst, h_dst);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: C=AB
 * @param {*} M A/C的行数
 * @param {*} N B/C的列数
 * @param {*} K A的列数，B的行数
 * @param {*} A 稠密矩阵
 * @param {*} B 稠密矩阵
 * @param {*} C 稠密矩阵（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Mult(const std::size_t M, const std::size_t N, const std::size_t K,
     const Kokkos::View<double **> &A, const Kokkos::View<double **> &B,
     Kokkos::View<double **> &C) {

  KokkosBlas::gemm("N", "N", static_cast<ScalarType>(1), A, B,
                   static_cast<ScalarType>(0), C);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: y=Ax
 * @param {*} M A的行数
 * @param {*} N A的列数
 * @param {*} A 稠密矩阵
 * @param {*} x 向量
 * @param {*} y 向量
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Mult(const std::size_t M, const std::size_t N, const Kokkos::View<double **> &A,
     const Kokkos::View<double *> &x, Kokkos::View<double *> &y) {

  KokkosBlas::gemv("N", static_cast<ScalarType>(1), A, x,
                   static_cast<ScalarType>(0), y);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: A = alpha*A
 * @param {*} alpha A的系数
 * @param {*} M A的行数
 * @param {*} N A的列数
 * @param {*} A 稠密矩阵（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Scal(const ScalarType alpha, const std::size_t M,
                                  const std::size_t N,
                                  Kokkos::View<double **> &A) {
  KokkosBlas::scal(A, alpha, A);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 获取矩阵元素A(i,j)
 * @param {*} i 行索引
 * @param {*} j 列索引
 * @param {*} M 行数
 * @param {*} N 列数
 * @param {*} A 稠密矩阵
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE ScalarType &
GetItem(const std::size_t i, const std::size_t j, const std::size_t M,
        const std::size_t N, Kokkos::View<double **> &A) {
  return A(i, j);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印矩阵A信息，一般用于调试
 * @param {*} M 行数
 * @param {*} N 列数
 * @param {*} A 稠密矩阵
 * @param {*} out 输出流（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Print(const std::size_t M, const std::size_t N,
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
    for (std::size_t j = 0; j < M - 1; ++j) {
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
