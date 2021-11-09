///
/// \file     scalar.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-01
/// \brief    标量用户接口，主要是为了衔接类似点积一类
///           操作的Device端实现。
///


#ifndef __CHIPSUM_NUMERIC_SCALAR_HPP__
#define __CHIPSUM_NUMERIC_SCALAR_HPP__

#include "../backend/backend.hpp"
#include "../chipsum_macro.h"

#include "numeric_traits.hpp"

#include "impl/scalar_serial_impl.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/scalar_kokkoskernels_impl.hpp"
// #include "impl/scalar_kokkos_impl.hpp"
#endif



namespace ChipSum {
namespace Numeric {


template <typename... Props>
class Scalar;


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



  CHIPSUM_DECLARED_FUNCTION
  ///
  /// \brief Scalar 构造函数
  ///
  Scalar() {
    ChipSum::Numeric::Impl::Scalar::create<ScalarType, SizeType>(__data);
  }


  CHIPSUM_DECLARED_FUNCTION
  ///
  /// \brief Scalar 构造函数
  ///
  Scalar(const ScalarType s) {
    ChipSum::Numeric::Impl::Scalar::create<ScalarType, SizeType>(s, __data);
  }


  CHIPSUM_FUNCTION_INLINE
  ///
  /// \brief GetData 获取后端底层，如Kokkos::View<double>
  /// \return 后端数据引用
  ///
  const_scalar_type_reference GetData() {
    return __data;
  }


  CHIPSUM_FUNCTION_INLINE

  ///
  /// \brief operator = 赋值操作符
  /// \param s POD类型
  /// \return Scalar数据对象引用（*this）
  ///
  Scalar& operator=(const ScalarType s) {
    ChipSum::Numeric::Impl::Scalar::deep_copy<ScalarType, SizeType>(s, __data);
    return *this;
  }

  CHIPSUM_FUNCTION_INLINE
  ///
  /// \brief operator = 赋值操作符
  /// \param s 后端数据类型
  /// \return
  ///
  Scalar& operator=(const scalar_type& s) {
    ChipSum::Numeric::Impl::Scalar::deep_copy<ScalarType, SizeType>(s, __data);
    return *this;
  }

  CHIPSUM_FUNCTION_INLINE
  ///
  /// \brief operator ()
  /// \return
  ///
  const ScalarType operator()() {
    return ChipSum::Numeric::Impl::Scalar::get_item<ScalarType, SizeType>(
        __data);
  }

  CHIPSUM_FUNCTION_INLINE
  ///
  /// \brief operator ScalarType
  ///
  operator ScalarType() const {
    return ChipSum::Numeric::Impl::Scalar::get_item<ScalarType, SizeType>(
        __data);
  }

  CHIPSUM_FUNCTION_INLINE
  ///
  /// \brief Print
  /// \param out
  ///
  void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::Scalar::print<ScalarType, SizeType>(__data, out);
  }
};

} // End namespace Numeric
} // End namespace ChipSum




typedef ChipSum::Numeric::Scalar<CSFloat, CSInt,
                                 ChipSum::Backend::DefaultBackend>
    Scalar;

#endif
