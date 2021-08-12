/*
 * @Description: 标量scalar的用户接口
 * @Version: 2.0
 * @Autor: lhl
 * @Date: 2021-08-09 12:27:29
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:44:04
 */


#ifndef __CHIPSUM_NUMERIC_SCALAR_HPP_
#define __CHIPSUM_NUMERIC_SCALAR_HPP_

#include "../backend/backend.hpp"
#include "../chipsum_macro.h"
#include "impl/scalar_kokkoskernels_impl.hpp"
#include "impl/scalar_serial_impl.hpp"
#include "numeric_traits.hpp"

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
  CHIPSUM_DECLARED_FUNCTION Scalar() {
    ChipSum::Numeric::Impl::Scalar::Create<ScalarType, SizeType>(__data);
  }

  CHIPSUM_DECLARED_FUNCTION Scalar(const ScalarType s) {
    ChipSum::Numeric::Impl::Scalar::Create<ScalarType, SizeType>(s, __data);
  }

  CHIPSUM_FUNCTION_INLINE const_scalar_type_reference GetData() {
    return __data;
  }

  CHIPSUM_FUNCTION_INLINE Scalar operator=(ScalarType &&s) {
    ChipSum::Numeric::Impl::Scalar::DeepCopy<ScalarType, SizeType>(s, __data);
  }

  CHIPSUM_FUNCTION_INLINE const ScalarType operator()() {
    return ChipSum::Numeric::Impl::Scalar::GetItem<ScalarType, SizeType>(
        __data);
  }

  CHIPSUM_FUNCTION_INLINE operator ScalarType() const {
    return ChipSum::Numeric::Impl::Scalar::GetItem<ScalarType, SizeType>(
        __data);
  }

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
