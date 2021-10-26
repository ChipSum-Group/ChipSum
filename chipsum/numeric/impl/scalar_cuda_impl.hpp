#ifndef __CHIPSUM_SCALAR_CUDA_IMPL_HPP__
#define __CHIPSUM_SCALAR_CUDA_IMPL_HPP__

#include <fstream>
#include <iostream>
#include <cassert>

#include <cuda_runtime.h>

#include "../../backend/backend.hpp"
#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

namespace ChipSum {
namespace Numeric {

template<typename ScalarType>
struct CudaScalar{
    ScalarType* val;


    ~CudaScalar(){
        cudaFree(val);
    }
};

template <typename ScalarType, typename SizeType, typename... Props>
struct Scalar_Traits<ScalarType, SizeType, ChipSum::Backend::Cuda, Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Cuda> {
  using scalar_type = CudaScalar<ScalarType>;

  using device_scalar_value_type = ScalarType;
};



namespace Impl {
namespace Scalar {
template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建标量
 * @param {View<ScalarType>} &r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Create(CudaScalar<ScalarType>& s) {

  cudaMalloc((void**)&s.val,sizeof (ScalarType));


}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建标量
 * @param {*} s 数据源
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Create(const ScalarType s, CudaScalar<ScalarType>& r) {
  cudaMalloc((void**)&r.val,sizeof (ScalarType));
  cudaMemcpy(r.val,&s,sizeof (ScalarType),cudaMemcpyHostToDevice);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 标量深拷贝
 * @param {*} s 数据源
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(const ScalarType s, CudaScalar<ScalarType>& r) {
  cudaMemcpy(r.val,&s,sizeof (ScalarType),cudaMemcpyHostToDevice);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 获取标量（by return）
 * @param {*} s 标量
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE ScalarType GetItem(CudaScalar<ScalarType>& r) {
  ScalarType s;
  cudaMemcpy(&s,r.val,sizeof (ScalarType),cudaMemcpyDeviceToHost);
  return s;
}




template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印标量，一般用于调试
 * @param {*} s 标量
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Print(const CudaScalar<ScalarType>& r, std::ostream &out) {

  ScalarType s;
  cudaMemcpy(&s,r.val,sizeof (ScalarType),cudaMemcpyDeviceToHost);
  cout << "scalar_serial: " << s << endl;
}
} // namespace Scalar
} // namespace Impl

} // namespace Numeric
} // namespace ChipSum


#endif // SCALAR_CUDA_IMPL_HPP
