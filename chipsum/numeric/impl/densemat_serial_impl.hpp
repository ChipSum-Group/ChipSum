/*
 * @Description: 稠密矩阵dense_matrix的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 15:57:04
 */

#ifndef __CHIPSUM_DENSEMAT_SERIAL_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_SERIAL_IMPL_HPP__

#include <cassert>
#include <fstream>
#include <vector>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"


/// 该头文件包含了稠密矩阵的算子实现。
/// 由于一些前期失误，dense matrix被设计成了类似blas的集合接口
/// 似乎是适合于处理小规模稠密矩阵运算的，这样的功能不太具有现实意义
/// 希望后续维护者能进行一些改进，改进建议如下：

/// 我希望dense matrix实际上依旧是面向代数系统求解的，
/// 所以它应当被设计成非常适合LU分解等操作的数据结构。
/// 也就是说矩阵的每一行要抽象为一个向量vector
/// 这样做的目的是为了能兼容如axpby操作，实现可分析扩展的算子组合
/// 例如高斯消元实际上就是不断地做axpby，最终形成三角阵
/// serial实现似乎不需要大规模变动底层数据结构，
/// 但kokkos后端的实现可能需要将view<ScalarType**>变为MultiVector，
/// 这部分工作我应该会亲自去实现，但后续的集成工作需要注意兼容。

/* 根据chipsum目前的工作来看，数据结构设计的重要性是远远大于算法创新的。 */

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
// 创建未初始化的矩阵
CHIPSUM_FUNCTION_INLINE void create(const std::size_t M, const std::size_t N,
                                    std::vector<ScalarType> &A) {

  A = std::vector<ScalarType>(M * N);
}

template <typename ScalarType, typename SizeType, typename... Props>
// 将POD数据填入矩阵
CHIPSUM_FUNCTION_INLINE void fill(const std::size_t M, const std::size_t N,
                                  ScalarType *src,
                                  std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(src, src + M * N);
}

template <typename ScalarType, typename SizeType, typename... Props>
// gemv
CHIPSUM_FUNCTION_INLINE void
mult(const std::size_t M, const std::size_t N, const std::vector<ScalarType> &A,
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
// gemm
CHIPSUM_FUNCTION_INLINE void
mult(const std::size_t M, const std::size_t N, const std::size_t K,
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
// A = alpha*A
CHIPSUM_FUNCTION_INLINE void scal(const ScalarType alpha, const std::size_t M,
                                  const std::size_t N,
                                  std::vector<ScalarType> &A) {
  assert(A.size() == M * N);
  for (std::size_t i = 0; i < A.size(); ++i) {
    A[i] *= alpha;
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
// 获取A(i,j)
CHIPSUM_FUNCTION_INLINE ScalarType &
get_item(const std::size_t i, const std::size_t j, const std::size_t M,
        const std::size_t N, std::vector<ScalarType> &A) {
  assert(i < M && j < N);

  return A[i * N + j];
}

template <typename ScalarType, typename SizeType, typename... Props>
// 打印矩阵信息
CHIPSUM_FUNCTION_INLINE void print(const std::size_t M, const std::size_t N,
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
