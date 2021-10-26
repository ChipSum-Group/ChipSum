/*
 * @Description: 标量scalar的用户接口
 * @Version: 2.0
 * @Autor: lhl
 * @Date: 2021-08-09 12:27:29
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-11 09:09:30
 */


#ifndef __CHIPSUM_NUMERIC_SCALAR_HPP_
#define __CHIPSUM_NUMERIC_SCALAR_HPP_

#include "../backend/backend.hpp"
#include "../chipsum_macro.h"

#include "numeric_traits.hpp"
#include "impl/scalar_cuda_impl.hpp"
#include "impl/scalar_serial_impl.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/scalar_kokkoskernels_impl.hpp"
// #include "impl/scalar_kokkos_impl.hpp"
#endif



namespace ChipSum {
namespace Numeric {

template <typename... Props> class Scalar;

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
class Scalar<ScalarType, SizeType, BackendType, Props...> {

public:
  using traits = Scalar_Traits<ScalarType, SizeType, BackendType, Props...>;

  using scalar_type = typename traits::scalar_type;

  using scalar_type_reference =
      typename std::add_lvalue_reference<scalar_type>::type;
  using const_scalar_type_reference =
      typename std::add_const<scalar_type_reference>::type;

  using device_scalar_type = typename traits::device_scalar_value_type;
  using device_scalar_type_reference =
      typename std::add_lvalue_reference<device_scalar_type>::type;
  ;

private:
  scalar_type __data;

public:

  /**
   * @description: 构造函数
   * @param {*}
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_DECLARED_FUNCTION Scalar() {
    ChipSum::Numeric::Impl::Scalar::Create<ScalarType, SizeType>(__data);
  }

  /**
   * @description: 构造函数（初始化）
   * @param {const ScalarType} s 初始值
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_DECLARED_FUNCTION Scalar(const ScalarType s) {
    ChipSum::Numeric::Impl::Scalar::Create<ScalarType, SizeType>(s, __data);
  }

  /**
   * @description: 获取数据
   * @param {*}
   * @return {const_scalar_type_reference} 数据
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE const_scalar_type_reference GetData() {
    return __data;
  }

  /**
   * @description: 赋值
   * @param {const ScalarType} 标量值
   * @return {Scalar} 
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE Scalar& operator=(const ScalarType s) {
    ChipSum::Numeric::Impl::Scalar::DeepCopy<ScalarType, SizeType>(s, __data);
    return *this;
  }

    /**
   * @description: 赋值
   * @param {scalar_type} 标量值
   * @return {Scalar} 
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE Scalar& operator=(const scalar_type& s) {
    ChipSum::Numeric::Impl::Scalar::DeepCopy<ScalarType, SizeType>(s, __data);
    return *this;
  }

  /**
   * @description: 获取标量值
   * @param {*} 
   * @return {const ScalarType} 标量值
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE const ScalarType operator()() {
    return ChipSum::Numeric::Impl::Scalar::GetItem<ScalarType, SizeType>(
        __data);
  }

  /**
   * @description: 隐式转换
   * @param {*} 
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE operator ScalarType() const {
    return ChipSum::Numeric::Impl::Scalar::GetItem<ScalarType, SizeType>(
        __data);
  }

  /**
   * @description: 打印（方便调试）
   * @param {ostream} &out 输出流
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::Scalar::Print<ScalarType, SizeType>(__data, out);
  }
};

} // End namespace Numeric
} // End namespace ChipSum




typedef ChipSum::Numeric::Scalar<double, std::size_t,
                                 ChipSum::Backend::DefaultBackend>
    Scalar;

#endif
