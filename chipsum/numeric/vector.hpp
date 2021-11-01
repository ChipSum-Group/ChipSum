///
/// \file     vector.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-10-27
///


/// \addtogroup 用户接口
#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <fstream>

#include "impl/vector_serial_impl.hpp"
#include "numeric_traits.hpp"
#include "scalar.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/vector_kokkoskernels_impl.hpp"
#include "impl/vector_dual_kokkos_impl.hpp"
#endif

namespace ChipSum {
namespace Numeric {


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


  CHIPSUM_DECLARED_FUNCTION
  Vector() = default;

  // 拷贝构造函数
  CHIPSUM_DECLARED_FUNCTION Vector(const Vector &) = default;

  // 拷贝构造函数
  CHIPSUM_DECLARED_FUNCTION Vector(Vector &&) = default;


  ///
  /// \brief Vector 构造函数
  /// \param size 向量维度
  ///
  CHIPSUM_DECLARED_FUNCTION Vector(const_size_type size) : __size(size) {
    ChipSum::Numeric::Impl::Vector::create<ScalarType, SizeType>(__size,
                                                                 __data);
  }


  ///
  /// \brief Vector 构造函数
  /// \param data 后端数据类型的源数据
  /// \param size 向量维度
  ///
  CHIPSUM_DECLARED_FUNCTION Vector(const vector_type &data,
                                   const_size_type size)
      : __size(size) {

    ChipSum::Numeric::Impl::Vector::create<ScalarType, SizeType>(__size,
                                                                 __data);

    ChipSum::Numeric::Impl::Vector::deep_copy<ScalarType, SizeType>(__data,
                                                                   data);
  }


  ///
  /// \brief Vector 构造函数
  /// \param data POD类型的数据源
  /// \param size 向量维度
  ///
  CHIPSUM_DECLARED_FUNCTION Vector(typename traits::nonconst_scalar_type *data,
                                   const_size_type size)
      : __size(size) {

    ChipSum::Numeric::Impl::Vector::create<ScalarType, SizeType>(data, size,
                                                                 __data);
  }


  ///
  /// \brief GetData 获取后端数据类型
  /// \return
  ///
  CHIPSUM_FUNCTION_INLINE const_vector_type_reference GetData() {
    return __data;
  }

  ///
  /// \brief GetSize 获取向量维度
  /// \return
  ///
  CHIPSUM_FUNCTION_INLINE const_size_type_reference GetSize() { return __size; }

  template <typename Arg>
  ///
  /// \brief Dot 莫版化向量点积
  /// \param v 右端项
  /// \param r 点积结果
  ///
  CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, Arg &r) {
    ChipSum::Numeric::Impl::Vector::dot<ScalarType, SizeType>(
        GetData(), v.GetData(), __size, r);
  }

  ///
  /// \brief Dot 向量点积
  /// \param v 右端项
  /// \param r 点积结果
  ///
  CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, scalar_type &r) {
    Dot(v, r.GetData());
  }

  template <typename RetType>
  ///
  /// \brief Dot 向量点积
  /// \param v 右端项
  /// \return 点积结果
  ///
  CHIPSUM_FUNCTION_INLINE RetType Dot(Vector &v) {
    RetType r;
    Dot(v, r);
    return r;
  }

  ///
  /// \brief Dot 向量点积
  /// \param v 右端项
  /// \return 点积结果
  ///
  CHIPSUM_FUNCTION_INLINE scalar_type Dot(Vector &v) {
    Scalar<ScalarType, SizeType, BackendType> r;
    Dot(v, r.GetData());
    return r;
  }

  ///
  /// \brief operator *  y = a*x
  /// \param a 系数
  /// \return y
  ///
  CHIPSUM_FUNCTION_INLINE Vector operator*(const ScalarType a) {

    Vector ret(__size);

    ChipSum::Numeric::Impl::Vector::scal<ScalarType, SizeType>(ret.GetData(), a,
                                                               GetData());
    return ret;
  }


  CHIPSUM_FUNCTION_INLINE Vector &operator=(const Vector &) = default;

  ///
  /// \brief operator * y = a*x
  /// \param a 系数(后端Scalar类型)
  /// \return y
  ///
  CHIPSUM_FUNCTION_INLINE Vector
  operator*(Scalar<ScalarType, SizeType, BackendType> &a) {

    Vector ret(__size);

    ChipSum::Numeric::Impl::Vector::scal<
        ScalarType, SizeType,
        typename Scalar<ScalarType, SizeType, BackendType>::device_scalar_type>(
        ret.GetData(), a.GetData(), __data);

    return ret;
  }

  ///
  /// \brief operator *=  x = a*x
  /// \param a 系数
  /// \return (*this)
  ///
  CHIPSUM_FUNCTION_INLINE Vector &operator*=(ScalarType a) {
    ChipSum::Numeric::Impl::Vector::scal<ScalarType, SizeType>(__data, a,
                                                               GetData());
    return *this;
  }

  ///
  /// \brief operator +  z=x+y
  /// \param y 右端项
  /// \return  结果
  ///
  CHIPSUM_FUNCTION_INLINE Vector operator+(Vector &y) {
    Vector ret(__data, __size);
    ChipSum::Numeric::Impl::Vector::axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(1), y.GetData(), ret.GetData());

    return ret;
  }

  ///
  /// \brief operator +=  x+=y
  /// \param y 右端项
  /// \return 结果
  ///
  CHIPSUM_FUNCTION_INLINE Vector &operator+=(Vector &y) {
    ChipSum::Numeric::Impl::Vector::axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(1), y.GetData(), __data);
    return *this;
  }

  ///
  /// \brief operator -  z=x-y
  /// \param y
  /// \return
  ///
  CHIPSUM_FUNCTION_INLINE Vector operator-(Vector &y) {

    Vector ret(__data, __size);
    ChipSum::Numeric::Impl::Vector::axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(-1), y.GetData(), ret.GetData());
    return ret;
  }

  ///
  /// \brief operator -  y=-x
  /// \return y
  ///
  CHIPSUM_FUNCTION_INLINE Vector operator-() {

    return static_cast<ScalarType>(-1) * (*this);
  }

  ///
  /// \brief operator -= x=x-y
  /// \param y
  /// \return
  ///
  CHIPSUM_FUNCTION_INLINE Vector &operator-=(Vector &y) {
    ChipSum::Numeric::Impl::Vector::axpy<ScalarType, SizeType>(
        static_cast<ScalarType>(-1), y.GetData(), __data);
    return *this;
  }

  ///
  /// \brief operator () x(i)
  /// \param i 下标索引
  /// \return
  ///
  CHIPSUM_FUNCTION_INLINE ScalarType &operator()(const SizeType i) {
    return ChipSum::Numeric::Impl::Vector::get_item<ScalarType, SizeType>(
        i, __data);
  }



  ///
  /// \brief Norm2 x的1范数
  /// \return 标量（POD）
  ///
  CHIPSUM_FUNCTION_INLINE ScalarType Norm1() {
    return ChipSum::Numeric::Impl::Vector::norm1<ScalarType, SizeType,
                                                 Props...>(__data);
  }


  ///
  /// \brief Norm2 x的2范数
  /// \return 标量（POD）
  ///
  CHIPSUM_FUNCTION_INLINE ScalarType Norm2() {
    return ChipSum::Numeric::Impl::Vector::norm2<ScalarType, SizeType,
                                                 Props...>(__data);
  }


  ///
  /// \brief NormInf x的Inf范数
  /// \return 标量（POD）
  ///
  CHIPSUM_FUNCTION_INLINE ScalarType NormInf() {
    return ChipSum::Numeric::Impl::Vector::norminf<ScalarType, SizeType,
                                                 Props...>(__data);
  }


  ///
  /// \brief Print 打印（一般用于调试）
  /// \param out 输出流
  ///
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::Vector::print<ScalarType, SizeType, Props...>(
        __data, out);
  }
};

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>


///
/// \brief operator * y=a*x
/// \param a 系数（POD）
/// \param x 向量
/// \return y
///
CHIPSUM_FUNCTION_INLINE Vector<ScalarType, SizeType, BackendType, Props...>
operator*(const ScalarType a,
          Vector<ScalarType, SizeType, BackendType, Props...> &x) {
  return x * a;
}

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
///
/// \brief operator * y=a*x
/// \param a 标量（后端数据类型）
/// \param x 向量x
/// \return
///
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
