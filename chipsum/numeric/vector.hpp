/*
 * @Description: 向量vector用户接口
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 14:32:38
 */

#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <fstream>

#include "impl/vector_kokkoskernels_impl.hpp"
#include "impl/vector_serial_impl.hpp"
#include "numeric_traits.hpp"
#include "scalar.hpp"

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

  CHIPSUM_DECLARED_FUNCTION Vector() = default;

  CHIPSUM_DECLARED_FUNCTION Vector(const Vector&) = default;

  CHIPSUM_DECLARED_FUNCTION Vector(Vector&&) = default;

  /**
   * @description:
   * @param {const_size_type} size
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION Vector(const_size_type size) : __size(size) {
    ChipSum::Numeric::Impl::Vector::Create<ScalarType, SizeType>(__size,
                                                                 __data);
  }

  /**
   * @description:
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
   * @description:
   * @param {typename} traits
   * @param {const_size_type} size
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION Vector(typename traits::nonconst_scalar_type *data,
                                   const_size_type size)
      : __size(size) {

    ChipSum::Numeric::Impl::Vector::Create<ScalarType, SizeType>(data, size,
                                                                 __data);
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE const_vector_type_reference GetData() {
    return __data;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE const_size_type_reference GetSize() { return __size; }

  template <typename RetType>
  /**
   * @description:
   * @param {Vector} &v
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE RetType Dot(Vector &v) {
    RetType r;
    ChipSum::Numeric::Impl::Vector::Dot<ScalarType, SizeType>(
        GetData(), v.GetData(), __size, r);
    return r;
  }

  /**
   * @description:
   * @param {Vector} &v
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE scalar_type Dot(Vector &v) {
    Scalar<ScalarType, SizeType, BackendType> r;
    Dot(v, r.GetData());
    return r;
  }

  template <typename Arg>
  /**
   * @description:
   * @param {Vector} &v
   * @param {Arg} &r
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, Arg &r) {
    ChipSum::Numeric::Impl::Vector::Dot<ScalarType, SizeType>(
        GetData(), v.GetData(), __size, r);
  }

  CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, scalar_type &r) {
    Dot(v, r.GetData());
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector operator*(const ScalarType s) {

    Vector ret(__size);
   
    ChipSum::Numeric::Impl::Vector::Scal<ScalarType, SizeType>(ret.GetData(), s,
                                                               GetData());


    return ret;
  }

  CHIPSUM_FUNCTION_INLINE Vector& operator=(const Vector&)=default;
  

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector
  operator*(Scalar<ScalarType, SizeType, BackendType>& s) {

    Vector ret(__size);

    ChipSum::Numeric::Impl::Vector::Scal<
        ScalarType, SizeType,
        typename Scalar<ScalarType, SizeType, BackendType>::device_scalar_type>(
        ret.GetData(), s.GetData(), __data);

    
    return ret;
  }


  
  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator*=(ScalarType s) {
    ChipSum::Numeric::Impl::Vector::Scal<ScalarType, SizeType>(__data, s,
                                                               GetData());
    return *this;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector operator+(Vector &x) {
    Vector ret(__data, __size);
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(1), x.GetData(), ret.GetData());

    return ret;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator+=(Vector &s) {
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(1), s.GetData(), __data);
    return *this;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector operator-(Vector &s) {

    Vector ret(__data, __size);
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(-1), s.GetData(), ret.GetData());
    return ret;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector operator-() {

    return static_cast<ScalarType>(-1) * (*this);
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE Vector &operator-=(Vector &s) {
    ChipSum::Numeric::Impl::Vector::Axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(-1), s.GetData(), __data);
    return *this;
  }

  /**
   * @description:
   * @param {const SizeType} i
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE ScalarType &operator()(const SizeType i) {
    return ChipSum::Numeric::Impl::Vector::GetItem<ScalarType, SizeType>(
        i, __data);
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE ScalarType Norm1() {
    return ChipSum::Numeric::Impl::Vector::Norm1<ScalarType, SizeType,
                                                 Props...>(__data);
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE ScalarType Norm2() {
    return ChipSum::Numeric::Impl::Vector::Norm2<ScalarType, SizeType,
                                                 Props...>(__data);
  }

  /**
   * @description:
   * @param {ostream} &out
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::Vector::Print<ScalarType, SizeType, Props...>(
        out, __data);
  }
};

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE Vector<ScalarType, SizeType, BackendType, Props...>
operator*(const ScalarType s,
          Vector<ScalarType, SizeType, BackendType, Props...> &v) {
  return v * s;
}

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE Vector<ScalarType, SizeType, BackendType, Props...>
operator*(Scalar<ScalarType, SizeType, BackendType, Props...> s,
          Vector<ScalarType, SizeType, BackendType, Props...> &v) {
  return v * s;
}


} // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::Vector<double, std::size_t,
                                 ChipSum::Backend::DefaultBackend>
    Vector;
typedef ChipSum::Numeric::Vector<double, std::size_t, ChipSum::Backend::Serial>
    SerialVector;

#endif // VECTOR_HPP
