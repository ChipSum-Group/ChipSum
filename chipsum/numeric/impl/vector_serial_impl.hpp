/*
 * @Description: 向量vector的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 16:11:28
 */

#ifndef __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__

#include <cassert>
#include <cmath>
#include <vector>

#include <fstream>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"


/// 很抱歉的是vector的设计出现了失误，而且出现在比较棘手的kokkos后端上
/// 这部分工作我自己不会去做，主要的原因是暂时没有时间。
/// 这部分工作将交给后来的参与者。

/// 对于vector的具体修改内容介绍我会放在vector_kokkoskernels_impl.hpp里面

namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Vector_Traits<ScalarType, SizeType, ChipSum::Backend::Serial, Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Serial> {
  using vector_type = typename std::vector<ScalarType>;
  using size_type = typename std::vector<ScalarType>::size_type;
  using device_scalar_value_type = ScalarType;
};

namespace Impl {

namespace Vector {

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void create(const SizeType n,
                                    std::vector<ScalarType> &dst) {
  dst.resize(n);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void create(const ScalarType *src, const SizeType n,
                                    std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(src, src + n);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void fill(const ScalarType val, const SizeType n,
                                  std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(n, val);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void dot(const std::vector<ScalarType> &x,

                                 const std::vector<ScalarType> &y,

                                 const SizeType &n, ScalarType &r) {

  assert(x.size() == y.size());

  for (SizeType i = 0; i < x.size(); ++i) {
    r += x[i] * y[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void scal(std::vector<ScalarType> &R,
                                  const ScalarType a,
                                  const std::vector<ScalarType> &X) {
  assert(R.size() == X.size());
  for (size_t i = 0; i < X.size(); ++i) {
    R[i] = a * X[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType norm1(const std::vector<ScalarType> &X) {
  ScalarType acc = 0.0;
  for (size_t i = 0; i < X.size(); ++i) {
    acc += X[i];
  }
  return acc;
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType norm2(const std::vector<ScalarType> &X) {
  ScalarType acc = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    acc += X[i] * X[i];
  }

  return std::sqrt(acc);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType norminf(const std::vector<ScalarType> &X) {
  ScalarType r = X[0];
  for (std::size_t i = 1; i < X.size(); ++i) {
    r = X[i]>r?X[i]:r;
  }
  return r;
}





template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void axpy(ScalarType a,
                                  const std::vector<ScalarType> &X,
                                  std::vector<ScalarType> &Y) {
  assert(X.size() == Y.size());
  for (std::size_t i = 0; i < Y.size(); ++i) {
    Y[i] += a * X[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void axpby(ScalarType a, std::vector<ScalarType> &X,
                                   ScalarType b, std::vector<ScalarType> &Y) {
  assert(X.size() == Y.size());
  for (std::size_t i = 0; i < Y.size(); ++i) {
    Y[i] = a * X[i] + b * Y[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void deep_copy(std::vector<ScalarType> &dst,
                                      const std::vector<ScalarType> &src) {
  dst.resize(src.size());
  for (std::size_t i = 0; i < dst.size(); ++i) {
    dst[i] = src[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void shallow_copy(std::vector<ScalarType> &dst,
                                         const std::vector<ScalarType> &src) {
  dst = src;
}


template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType &get_item(const std::size_t index,
                                            std::vector<ScalarType> &vec) {

  assert(index < vec.size());
  return vec[index];
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void print(const std::vector<ScalarType> &vec,
                                   std::ostream &out) {

  out << " [";
  for (std::size_t i = 0; i < vec.size() - 1; ++i) {
    out << vec[i] << ", ";
  }

  out << vec[vec.size() - 1] << "]" << std::endl;
}


} // End namespace Vector
} // End namespace Impl

} // End namespace Numeric

} // End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
