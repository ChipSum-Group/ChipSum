/*
 * @Description: 稠密矩阵dense_matrix的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-19 09:10:23
 */

#ifndef __CHIPSUM_DENSEMAT_SERIAL_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_SERIAL_IMPL_HPP__

#include <cassert>
#include <fstream>
#include <vector>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

namespace ChipSum {

namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct DenseMatrix_Traits<ScalarType, SizeType, ChipSum::Backend::Serial,
                          Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Serial,
                             Props...> {

  using matrix_type = std::vector<ScalarType>;

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
                                    std::vector<ScalarType> &A) {

  A = std::vector<ScalarType>(M * N);
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
                                  std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(src, src + M * N);
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
Mult(const std::size_t M, const std::size_t N, const std::vector<ScalarType> &A,
     const std::vector<ScalarType> &x, std::vector<ScalarType> &b) {

  for (std::size_t i = 0; i < M; ++i)
    b[i] = 0;

  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      b[i] += A[i * N + j] * x[j];
    }
  }
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
     const std::vector<ScalarType> &A, const std::vector<ScalarType> &B,
     std::vector<ScalarType> &C) {

  

  for (std::size_t i = 0; i < C.size(); ++i)
    C[i] = 0;

  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < K; ++j) {
      ScalarType Aik = A[i * K + j];
      for (std::size_t k = 0; k < N; ++k) {
        C[i * N + k] += Aik * B[j * K + k];
      }
    }
  }
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
                                  std::vector<ScalarType> &A) {
  assert(A.size() == M * N);
  for (std::size_t i = 0; i < A.size(); ++i) {
    A[i] *= alpha;
  }
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
        const std::size_t N, std::vector<ScalarType> &A) {
  assert(i < M && j < N);

  return A[i * N + j];
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印矩阵A信息，一般用于调试
 * @param {*} M 行数
 * @param {*} N 列数
 * @param {*} A 稠密矩阵
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Print(const std::size_t M, const std::size_t N,
                                   std::vector<ScalarType> &A,
                                   std::ostream &out) {

  cout << "dense_mat_serial("
  << M <<","
  << N <<")"
  << ":" << endl;

  for (std::size_t i = 0; i < M; ++i) {

    out << " "
        << "[";
    for (std::size_t j = 0; j < N - 1; ++j) {
      out << A[i * N + j] << ", ";
    }
    out << A[i * N + M - 1] << "]" << endl;
  }
  out << endl;
}

} // namespace DenseMat
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
