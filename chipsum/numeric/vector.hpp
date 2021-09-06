/*
 * @Description: 向量vector用户接口
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-23 09:30:24
 */

#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <fstream>

#include "impl/vector_serial_impl.hpp"
#include "numeric_traits.hpp"
#include "scalar.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/vector_kokkoskernels_impl.hpp"
#endif

namespace ChipSum {
namespace Numeric {

template <typename... Props> class Vector;

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
class Vector<ScalarType, SizeType, BackendType, Props...> {

public:
  using traits = Vector_Traits<ScalarType, SizeType, BackendType, Props...>;

  using vector_type = typename traits::vector_type;

  using size_type = typename traits::size_type;
  using const_size_type = typename traits::const_size_type;
  using size_type_reference =
      typename std::add_lvalue_reference<size_type>::type;
  using const_size_type_reference =
      typename std::add_const<size_type_reference>::type;

  using vector_type_reference =
      typename std::add_lvalue_reference<vector_type>::type;
  using const_vector_type_reference =
      typename std::add_const<vector_type_reference>::type;

  using device_scalar_type = typename traits::device_scalar_value_type;
  using device_scalar_type_reference =
      typename std::add_lvalue_reference<device_scalar_type>::type;

  using scalar_type = Scalar<ScalarType, SizeType, BackendType>;

private:
  vector_type __data;
  size_type __size;

private:
  //    friend class DenseMatrix<ScalarType,SizeType,BackendType,Props...>;

  //    void privateSample(){}

protected:
  //    void _protectedSample(){}

public:
  //    void PublicSample(){}

  /**
   * @description: 构造函数
   * @param {*}
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_DECLARED_FUNCTION Vector() = default;

  // 拷贝构造函数
  CHIPSUM_DECLARED_FUNCTION Vector(const Vector &) = default;

  // 拷贝构造函数
  CHIPSUM_DECLARED_FUNCTION Vector(Vector &&) = default;

  /**
   * @description: 创建长度为size的Vector
   * @param {const_size_type} size
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION Vector(const_size_type size) : __size(size) {
    ChipSum::Numeric::Impl::Vector::Create<ScalarType, SizeType>(__size,
                                                                 __data);
  }

  /**
   * @description: 构造函数
   * @param {const vector_type} &data
   * @param {const_size_type} size
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION Vector(const vector_type &data,
                                   const_size_type size)
      : __size(size) {

    ChipSum::Numeric::Impl::Vector::Create<ScalarType, SizeType>(__size,
                                                                 __data);

    ChipSum::Numeric::Impl::Vector::DeepCopy<ScalarType, SizeType>(__data,
                                                                   data);
  }

  /**
   * @description: 拷贝构造函数
   * @param {nonconst_scalar_type} data POD数据
   * @param {const_size_type} size 向量长度
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION Vector(typename traits::nonconst_scalar_type *data,
                                   const_size_type size)
      : __size(size) {

    ChipSum::Numeric::Impl::Vector::Create<ScalarType, SizeType>(data, size,
                                                                 __data);
  }

  /**
   * @description: 获取向量数据
   * @param {*}
   * @return {const_vector_type_reference} 数据
   */
  CHIPSUM_FUNCTION_INLINE const_vector_type_reference GetData() {
    return __data;
  }

  /**
   * @description: 获取向量长度
   * @param {*}
   * @return {const_size_type_reference} 向量长度
   */
  CHIPSUM_FUNCTION_INLINE const_size_type_reference GetSize() { return __size; }

  template <typename Arg>
  /**
   * @description: dot运算
   * @param {Vector} v 向量
   * @param {Arg} r 标量
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, Arg &r) {
    ChipSum::Numeric::Impl::Vector::Dot<ScalarType, SizeType>(
        GetData(), v.GetData(), __size, r);
  }

  /**
   * @description: dot运算
   * @param {Vector} v 向量
   * @param {scalar_type} r 标量
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, scalar_type &r) {
    Dot(v, r.GetData());
  }

  template <typename RetType>
  /**
   * @description: dot运算
   * @param {Vector} v 向量
   * @return {RetType} 标量
   */
  CHIPSUM_FUNCTION_INLINE RetType Dot(Vector &v) {
    RetType r;
    Dot(v, r);
    return r;
  }

  /**
   * @description: dot运算
   * @param {Vector} v 向量
   * @return {scalar_type} 标量
   */
  CHIPSUM_FUNCTION_INLINE scalar_type Dot(Vector &v) {
    Scalar<ScalarType, SizeType, BackendType> r;
    Dot(v, r.GetData());
    return r;
  }

  /**
   * @description: y=a*x
   * @param {const ScalarType} a 系数
   * @return {Vector} 结果
   */
  CHIPSUM_FUNCTION_INLINE Vector operator*(const ScalarType a) {

    Vector ret(__size);

    ChipSum::Numeric::Impl::Vector::Scal<ScalarType, SizeType>(ret.GetData(), a,
                                                               GetData());
    return ret;
  }

  /**
   * @description: 赋值拷贝运算符
   * @param {const Vector&}
   * @return {Vector&}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator=(const Vector &) = default;

  /**
   * @description: y=a*x
   * @param {Scalar<...>} a 系数
   * @return {Vector} 向量
   */
  CHIPSUM_FUNCTION_INLINE Vector
  operator*(Scalar<ScalarType, SizeType, BackendType> &a) {

    Vector ret(__size);

    ChipSum::Numeric::Impl::Vector::Scal<
        ScalarType, SizeType,
        typename Scalar<ScalarType, SizeType, BackendType>::device_scalar_type>(
        ret.GetData(), a.GetData(), __data);

    return ret;
  }

  /**
   * @description: x*=a
   * @param {ScalarType} a 系数
   * @return {Vector&} 向量
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator*=(ScalarType a) {
    ChipSum::Numeric::Impl::Vector::Scal<ScalarType, SizeType>(__data, a,
                                                               GetData());
    return *this;
  }

  /**
   * @description: c=a+b
   * @param {Vector&} b
   * @return {Vector} 向量
   */
  CHIPSUM_FUNCTION_INLINE Vector operator+(Vector &b) {
    Vector ret(__data, __size);
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(1), b.GetData(), ret.GetData());

    return ret;
  }

  /**
   * @description: x+=y
   * @param {Vector&} y 向量
   * @return {Vector&} 向量
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator+=(Vector &y) {
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(1), y.GetData(), __data);
    return *this;
  }

  /**
   * @description: c = a-b
   * @param {Vector&} b 向量
   * @return {Vector} 向量
   */
  CHIPSUM_FUNCTION_INLINE Vector operator-(Vector &b) {

    Vector ret(__data, __size);
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(-1), b.GetData(), ret.GetData());
    return ret;
  }

  /**
   * @description: -x
   * @param {*}
   * @return {Vector} 结果
   */
  CHIPSUM_FUNCTION_INLINE Vector operator-() {

    return static_cast<ScalarType>(-1) * (*this);
  }

  /**
   * @description: x-=y
   * @param {Vector} y 向量
   * @return {Vector} 向量
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator-=(Vector &y) {
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(-1), y.GetData(), __data);
    return *this;
  }

  /**
   * @description: 下标取值
   * @param {const SizeType} i 下标索引
   * @return {ScalarType &} 标量（POD）
   */
  CHIPSUM_FUNCTION_INLINE ScalarType &operator()(const SizeType i) {
    return ChipSum::Numeric::Impl::Vector::GetItem<ScalarType, SizeType>(
        i, __data);
  }

  /**
   * @description: x的1范数
   * @param {device_scalar_type_reference} r 1范数结果
   * @return {*} 标量（POD）
   */
  CHIPSUM_FUNCTION_INLINE void Norm1(device_scalar_type_reference r) {
    ChipSum::Numeric::Impl::Vector::Norm1<ScalarType, SizeType, Props...>(
        r, __data);
  }

  /**
   * @description: x的2范数
   * @param {device_scalar_type_reference} r 2范数结果
   * @return {*} 标量（POD）
   */
  CHIPSUM_FUNCTION_INLINE void Norm2(device_scalar_type_reference r) {
    ChipSum::Numeric::Impl::Vector::Norm2<ScalarType, SizeType, Props...>(
        r, __data);
  }

  /**
   * @description: x的1范数
   * @param {*}
   * @return {*} 标量（POD）
   */
  CHIPSUM_FUNCTION_INLINE ScalarType Norm1() {
    return ChipSum::Numeric::Impl::Vector::Norm1<ScalarType, SizeType,
                                                 Props...>(__data);
  }

  /**
   * @description: x的2范数
   * @param {*}
   * @return {*} 标量（POD）
   */
  CHIPSUM_FUNCTION_INLINE ScalarType Norm2() {
    return ChipSum::Numeric::Impl::Vector::Norm2<ScalarType, SizeType,
                                                 Props...>(__data);
  }

  /**
   * @description: 打印（一般用于调试）
   * @param {std::ostream&} out 输出流
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::Vector::Print<ScalarType, SizeType, Props...>(
        __data, out);
  }
};

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
/**
 * @description: y=a*x
 * @param {const ScalarType} 系数（POD）
 * @param {Vector} 向量
 * @return {*} 向量
 */
CHIPSUM_FUNCTION_INLINE Vector<ScalarType, SizeType, BackendType, Props...>
operator*(const ScalarType a,
          Vector<ScalarType, SizeType, BackendType, Props...> &x) {
  return x * a;
}

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
/**
 * @description: y=a*x
 * @param {Scalar} 标量（ChipSum）
 * @param {Vector} 向量
 * @return {*} 向量
 */
CHIPSUM_FUNCTION_INLINE Vector<ScalarType, SizeType, BackendType, Props...>
operator*(Scalar<ScalarType, SizeType, BackendType, Props...> &a,
          Vector<ScalarType, SizeType, BackendType, Props...> &x) {
  return x * a;
}

} // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::Vector<double, std::size_t,
                                 ChipSum::Backend::DefaultBackend>
    Vector;
typedef ChipSum::Numeric::Vector<double, std::size_t, ChipSum::Backend::Serial>
    SerialVector;

#endif // VECTOR_HPP
