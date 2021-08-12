/*
 * @Description: 稠密矩阵dense_matrix的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:42:58
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
 * @brief Create
 * @param M
 * @param N
 * @param mat
 */
CHIPSUM_FUNCTION_INLINE void Create(const std::size_t M, const std::size_t N,
                                    std::vector<ScalarType> &mat) {

  mat = std::vector<ScalarType>(M * N);
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
                                  std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(src, src + M * N);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Mult
 * @param M
 * @param N
 * @param A
 * @param x
 * @param b
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
 * @brief Scal
 * @param alpha
 * @param M
 * @param N
 * @param mat
 */
CHIPSUM_FUNCTION_INLINE void Scal(ScalarType alpha, const std::size_t M,
                                  const std::size_t N,
                                  std::vector<ScalarType> &mat) {
  assert(mat.size() == M * N);
  for (std::size_t i = 0; i < mat.size(); ++i) {
    mat[i] *= alpha;
  }
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
        const std::size_t N, std::vector<ScalarType> &mat) {
  assert(i < M && j < N);

  return mat[i * N + j];
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
                                   std::vector<ScalarType> &mat,
                                   std::ostream &out) {

  for (std::size_t i = 0; i < M; ++i) {

    out << " "
        << "[";
    for (std::size_t j = 0; j < M - 1; ++j) {
      out << mat[i * N + j] << ", ";
    }
    out << mat[i * N + M - 1] << "]" << endl;
  }
  out << endl;
}

} // namespace DenseMat
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
