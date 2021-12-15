///
/// \file     vector.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-01
/// \brief    向量类用户接口


#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <fstream>

#include "numeric_traits.hpp"
#include "scalar.hpp"



#include "../common/data_types.hpp"

#include "impl/vector_impl.hpp"

#include "impl/impl_abstract.hpp"

namespace ChipSum {
namespace Numeric {




template <typename ...Props>
class Vector {

public:
    using traits = Vector_Traits<Props...>;

    using vector_type = typename traits::vector_type;

    using size_type = typename traits::size_type;
//    using const_size_type = typename traits::const_size_type;
    using size_type_ref =
    typename ::std::add_lvalue_reference<size_type>::type;
    using const_size_type_ref =
    typename ::std::add_const<size_type_ref>::type;

    using vector_type_ref =
    typename ::std::add_lvalue_reference<vector_type>::type;
    using const_vector_type_ref =
    typename ::std::add_const<vector_type_ref>::type;

    using scalar_type = Scalar<Props...>;


    using value_type = typename traits::value_type;
    using value_type_ref =
    typename ::std::add_lvalue_reference<value_type>::type;
    using const_value_type_ref =
    typename ::std::add_const<value_type_ref>::type;

private:
    vector_type __data;
    size_type __size;



protected:


public:



    CHIPSUM_DECLARED_FUNCTION
    Vector() = default;

    // 拷贝构造函数
    CHIPSUM_DECLARED_FUNCTION Vector(const Vector &) = default;

    // 拷贝构造函数
    CHIPSUM_DECLARED_FUNCTION Vector(Vector &&) = default;






    template<typename ST,typename ...Args>

    CHIPSUM_DECLARED_FUNCTION Vector(const ST& size,Args... args)
        :__size(size)
    {
        ChipSum::Numeric::Impl::Vector::create(__data,size,args...);
    }





    CHIPSUM_FUNCTION_INLINE void DeepCopy(Vector& y){
        ChipSum::Numeric::Impl::Vector::deep_copy(__data,y.GetData());
    }


    ///
    /// \brief GetData 获取后端数据类型
    /// \return
    ///
    CHIPSUM_FUNCTION_INLINE vector_type& GetData() {
        return __data;
    }

    ///
    /// \brief GetSize 获取向量维度
    /// \return
    ///
    CHIPSUM_FUNCTION_INLINE const_size_type_ref GetSize() { return __size; }

    template <typename Arg>

    CHIPSUM_FUNCTION_INLINE void Dot(Vector &v, Arg& r) {
        ChipSum::Numeric::Impl::Vector::dot(
                    GetData(), v.GetData(),r);
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
    CHIPSUM_FUNCTION_INLINE value_type Dot(Vector &v) {

        return ChipSum::Numeric::Impl::Vector::dot(__data,v.GetData());
    }

    ///
    /// \brief Dot 向量点积
    /// \param v 右端项
    /// \return 点积结果
    ///
    CHIPSUM_FUNCTION_INLINE scalar_type Dot(Vector &v) {
        Scalar<Props...> r;
        Dot(v, r.GetData());
        return r;
    }

    ///
    /// \brief operator *  y = a*x
    /// \param a 系数
    /// \return y
    ///
    CHIPSUM_FUNCTION_INLINE Vector operator*(value_type a) {

        Vector ret(__size);

        ChipSum::Numeric::Impl::Vector::scal(__data, a, ret.GetData());
        return ret;
    }



    template<typename VT>
    ///
    /// \brief Mult y=a*(*this)
    /// \param a
    /// \param y
    ///
    CHIPSUM_FUNCTION_INLINE void Mult(const VT& a,Vector& y) {

        ChipSum::Numeric::Impl::Vector::scal(__data, a, y.GetData());
    }


    CHIPSUM_FUNCTION_INLINE Vector &operator=(Vector& y){

        DeepCopy(y);
        return (*this);

    }









    ///
    /// \brief operator * y = a*x
    /// \param a 系数(后端Scalar类型)
    /// \return y
    ///
    CHIPSUM_FUNCTION_INLINE void
    Mult(Scalar<Props...> &a,Vector& y) {



        ChipSum::Numeric::Impl::Vector::scal(
                     __data,a.GetData(),y.GetData());


    }

    ///
    /// \brief operator * y = a*x
    /// \param a 系数(后端Scalar类型)
    /// \return y
    ///
    CHIPSUM_FUNCTION_INLINE Vector&
    operator*=(Scalar<Props...> &a) {
        ChipSum::Numeric::Impl::Vector::scal(
                     __data,a.GetData(),__data.GetData());

        return (*this);
    }

    ///
    /// \brief operator *=  x = a*x
    /// \param a 系数
    /// \return (*this)
    ///
    CHIPSUM_FUNCTION_INLINE Vector &operator*=(value_type a) {
        ChipSum::Numeric::Impl::Vector::scal(__data, a,
                                             __data);
        return *this;
    }



    ///
    /// \brief operator +  z=x+y
    /// \param y 右端项
    /// \return  结果
    ///
    CHIPSUM_FUNCTION_INLINE Vector Add(Vector &y) {

        ChipSum::Numeric::Impl::Vector::add(
                    __data, y.GetData());


    }






    ///
    /// \brief operator +=  x+=y
    /// \param y 右端项
    /// \return 结果
    ///
    CHIPSUM_FUNCTION_INLINE Vector &operator+=(Vector &y) {
        ChipSum::Numeric::Impl::Vector::add(
                    y.GetData(), __data);
        return *this;
    }



    ///
    /// \brief operator -  z=x-y
    /// \param y
    /// \return
    ///
    CHIPSUM_FUNCTION_INLINE void Sub(Vector &y) {



        ChipSum::Numeric::Impl::Vector::axpby(
                    y.GetData(), GetData(),-1,1);

    }

    ///
    /// \brief Sub 减法接口 z = (*this)-y
    /// \param y 输入 y
    /// \param z 输出 z
    ///
    CHIPSUM_FUNCTION_INLINE void Sub(Vector &y,Vector &z) {




        ChipSum::Numeric::Impl::Vector::deep_copy(z.GetData(),__data);


        ChipSum::Numeric::Impl::Vector::axpby(
                    y.GetData(), z.GetData(),-1,1);

    }

    ///
    /// \brief operator -  y=-x
    /// \return y
    ///
    CHIPSUM_FUNCTION_INLINE Vector operator-() {

        return -1 * (*this);
    }

    ///
    /// \brief operator -  y=-x
    /// \return y
    ///
    CHIPSUM_FUNCTION_INLINE Vector Neg() {

        return Mult(-1);
    }

    ///
    /// \brief operator -= x=x-y
    /// \param y
    /// \return
    ///
    CHIPSUM_FUNCTION_INLINE Vector &operator-=(Vector &y) {
        ChipSum::Numeric::Impl::Vector::axpby(
                    -1, y.GetData(), __data);
        return *this;
    }

    ///
    /// \brief operator () x(i)
    /// \param i 下标索引
    /// \return
    ///
    CHIPSUM_FUNCTION_INLINE value_type_ref operator()(size_type i) {
        return ChipSum::Numeric::Impl::Vector::get_item(
                    i, __data);
    }



    ///
    /// \brief Norm2 x的1范数
    /// \return 标量（POD）
    ///
    CHIPSUM_FUNCTION_INLINE auto Norm1() {
        return ChipSum::Numeric::Impl::Vector::norm1(__data);
    }


    template<typename RefType>
    ///
    /// \brief Norm1 1范数
    /// \param r [OUT]
    ///
    CHIPSUM_FUNCTION_INLINE void Norm1(RefType& r) {
        ChipSum::Numeric::Impl::Vector::norm1(__data,r);
    }

    CHIPSUM_FUNCTION_INLINE void Norm1(scalar_type& r) {
        ChipSum::Numeric::Impl::Vector::norm1(__data,r.GetData());
    }


    ///
    /// \brief Norm2 x的2范数
    /// \return 标量（POD）
    ///
    CHIPSUM_FUNCTION_INLINE value_type Norm2() {
        return ChipSum::Numeric::Impl::Vector::norm2(__data);
    }


    template<typename RefType>
    ///
    /// \brief Norm1 2范数
    /// \param r [OUT]
    ///
    CHIPSUM_FUNCTION_INLINE void Norm2(RefType& r) {
        ChipSum::Numeric::Impl::Vector::norm2(__data,r);
    }

    CHIPSUM_FUNCTION_INLINE void Norm2(scalar_type& r) {
        ChipSum::Numeric::Impl::Vector::norm2(__data,r.GetData());
    }
    ///
    /// \brief NormInf x的Inf范数
    /// \return 标量（POD）
    ///
    CHIPSUM_FUNCTION_INLINE value_type NormInf() {
        return ChipSum::Numeric::Impl::Vector::norminf(__data);
    }

    template<typename RefType>
    ///
    /// \brief Norm1 2范数
    /// \param r [OUT]
    ///
    CHIPSUM_FUNCTION_INLINE void NormInf(RefType& r) {
        ChipSum::Numeric::Impl::Vector::norminf(__data,r);
    }

    CHIPSUM_FUNCTION_INLINE void NormInf(scalar_type& r) {
        ChipSum::Numeric::Impl::Vector::norminf(__data,r.GetData());
    }

    template<typename OStreamT = std::ostream>
    ///
    /// \brief Print 打印（一般用于调试）
    /// \param out 输出流
    ///
    CHIPSUM_FUNCTION_INLINE void Print(OStreamT &out = std::cout) {
        ChipSum::Numeric::Impl::Vector::print(
                    __data, out);
    }


    template<typename ...Args>
    ///
    /// \brief AXPBY 模板化接口
    /// \param y  输出，可能为 y += x  y += a*x y = a*x+b*y
    /// \param args
    ///
    CHIPSUM_FUNCTION_INLINE void AXPBY(Vector& y,Args ...args)
    {

        ChipSum::Numeric::Impl::Vector::axpby(GetData(),y.GetData(),args...);
    }
};

template <typename ValueType,typename... Props>
///
/// \brief operator * y=a*x
/// \param a 系数（POD）
/// \param x 向量
/// \return y
///
CHIPSUM_FUNCTION_INLINE Vector< ValueType,Props...>
operator*(const ValueType& a,
          Vector<ValueType,Props...> &x) {
    return x * a;
}




template <typename... Props>
///
/// \brief operator * y=a*x
/// \param a 标量（后端数据类型）
/// \param x 向量x
/// \return
///
CHIPSUM_FUNCTION_INLINE Vector<Props...>
operator*(Scalar<Props...> &a,
          Vector<Props...> &x) {
    return x * a;
}

} // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::Vector<CSFloat,
ChipSum::Backend::DefaultBackend>
Vector;

typedef ChipSum::Numeric::Vector<CSFloat,
ChipSum::Backend::Serial>
SerialVector;

#endif // VECTOR_HPP
