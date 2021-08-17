/*
 * @Description: 向量vector的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-16 15:40:31
 */

#ifndef __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__

#include <cassert>
#include <cmath>
#include <vector>

#include <fstream>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

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
/**
 * @description: 创建未初始化的向量
 * @param {*} n 向量长度 
 * @param {*} x 向量（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,
                                    std::vector<ScalarType> &dst) {
  dst.resize(n);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建初始化的向量
 * @param {*} src POD数据源
 * @param {*} n 向量长度
 * @param {*} x 向量（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Create(const ScalarType *src, const SizeType n,
                                    std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(src, src + n);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Fill：利用POD数据进行向量填充
 * @param src：POD数据源
 * @param n：向量维度
 * @param dst：目标向量
 */
CHIPSUM_FUNCTION_INLINE void Fill(const ScalarType val, const SizeType n,
                                  std::vector<ScalarType> &dst) {
  dst = std::vector<ScalarType>(n, val);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 向量dot运算（by reference）
 * @param {*} a 向量
 * @param {*} b 向量
 * @param {*} n 向量长度 
 * @param {*} r 结果（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Dot(const std::vector<ScalarType> &x,

                                 const std::vector<ScalarType> &y,

                                 const SizeType &n, ScalarType &r) {

  assert(x.size() == y.size());

  for (SizeType i = 0; i < x.size(); ++i) {
    r += x[i] * y[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: R=a*X
 * @param {*} R Scal结果（out）
 * @param {*} a X的系数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Scal(std::vector<ScalarType> &R,
                                  const ScalarType a,
                                  const std::vector<ScalarType> &X) {
  assert(R.size() == X.size());
  for (size_t i = 0; i < X.size(); ++i) {
    R[i] = a * X[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: X的1范数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm1(const std::vector<ScalarType> &X) {
  ScalarType acc = 0.0;
  for (size_t i = 0; i < X.size(); ++i) {
    acc += X[i];
  }
  return acc;
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: X的2范数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm2(const std::vector<ScalarType> &X) {
  ScalarType acc = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i) {
    acc += X[i] * X[i];
  }

  return std::sqrt(acc);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: Y=Y+a*X
 * @param {*} a 系数
 * @param {*} X 向量
 * @param {*} Y 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Axpy(ScalarType a,
                                  const std::vector<ScalarType> &X,
                                  std::vector<ScalarType> &Y) {
  assert(X.size() == Y.size());
  for (std::size_t i = 0; i < Y.size(); ++i) {
    Y[i] += a * X[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: Y=b*Y+a*X
 * @param {*} a 系数
 * @param {*} X 向量
 * @param {*} b 系数
 * @param {*} Y 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Axpby(ScalarType a, std::vector<ScalarType> &X,
                                   ScalarType b, std::vector<ScalarType> &Y) {
  assert(X.size() == Y.size());
  for (std::size_t i = 0; i < Y.size(); ++i) {
    Y[i] = a * X[i] + b * Y[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 深拷贝
 * @param {*} dst 目标数据
 * @param {*} src 原数据
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(std::vector<ScalarType> &dst,
                                      const std::vector<ScalarType> &src) {
  dst.resize(src.size());
  for (std::size_t i = 0; i < dst.size(); ++i) {
    dst[i] = src[i];
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 浅拷贝
 * @param {*} dst 目标数据
 * @param {*} src 原数据
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void ShallowCopy(std::vector<ScalarType> &dst,
                                         const std::vector<ScalarType> &src) {
  dst = src;
}


template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 获取向量元素
 * @param {*} index 索引
 * @param {*} vec 向量
 * @return {*} 元素值
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE ScalarType &GetItem(const std::size_t index,
                                            std::vector<ScalarType> &vec) {

  assert(index < vec.size());
  return vec[index];
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印向量，一般用于调试
 * @param {*} vec 向量
 * @param {*} out 输出流
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Print(const std::vector<ScalarType> &vec,
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
