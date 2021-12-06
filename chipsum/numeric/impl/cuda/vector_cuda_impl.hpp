/*
 * @Description: 向量vector的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 16:16:28
 */

#ifndef __CHIPSUM_VECTOR_CUDA_IMPL_HPP__
#define __CHIPSUM_VECTOR_CUDA_IMPL_HPP__

#include <cassert>
#include <vector>


#include <fstream>
#include <cuda_runtime.h>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

namespace ChipSum {
namespace Numeric {


template<typename ScalarType>
struct CudaVector{
  ScalarType* d_val;
  std::vector<ScalarType> h_val;

  std::size_t size;

  ~CudaVector(){
    cudaFree(d_val);  
  }
};

template <typename ScalarType, typename SizeType, typename... Props>
struct Vector_Traits<ScalarType, SizeType, ChipSum::Backend::Cuda, Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Cuda> {
  using vector_type = CudaVector;
  using size_type = typename std::vector<ScalarType>::size_type;
  using device_scalar_value_type = ScalarType;
};

namespace Impl {

namespace Vector {

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void create(const SizeType n,
                                    CudaVector<ScalarType> &dst) {
  
  std::size_t vec_size = static_cast<std::size_t>(n);
  cudaMalloc((void**)&dst.d_val,vec_size*sizeof(ScalarType));
  dst.size = vec_size;
  dst.h_val.resize(vec_size,0);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void create(const ScalarType *src, const SizeType n,
                                    CudaVector<ScalarType> &dst) {
  std::size_t vec_size = static_cast<std::size_t>(n);
  dst.h_val = std::vector<ScalarType>(src, src + vec_size);
  cudaMalloc((void**)&dst.d_val,vec_size*sizeof(ScalarType));
  cudaMemcpy(dst.d_val,src,vec_size*sizeof(ScalarType),cudaMemcpyHostToDevice);

}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void fill(const ScalarType val, const SizeType n,
                                  CudaVector<ScalarType> &dst) {
  std::size_t vec_size = static_cast<std::size_t>(n);
  dst.h_val = std::vector<ScalarType>(vec_size,val);
  cudaMalloc((void**)&dst.d_val,vec_size*sizeof(ScalarType));
  cudaMemcpy(dst.d_val,dst.h_val.data(),vec_size*sizeof(ScalarType),cudaMemcpyHostToDevice);
}


template <typename ScalarType>
__global__ void dot(ScalarType* x,ScalarType* y,std::size_t n){
  
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
