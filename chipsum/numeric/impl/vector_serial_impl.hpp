/* * * * * * * * * * * * * * * * * * * * *
*   File:     vector_serial_impl.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__

#include <vector>
#include <cmath>
#include <cassert>

#include <fstream>



#include "../numeric_traits.hpp"
#include "../../chipsum_macro.h"


namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,ChipSum::Backend::CPUSerialBackend,Props...>:
        public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::CPUSerialBackend>
{
    using vector_type = typename std::vector<ScalarType>;
    using size_type = typename std::vector<ScalarType>::size_type;
    using device_scalar_value_type = ScalarType;

};




namespace Impl {

namespace Vector {

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Create：创建向量数据，主要是申请内存空间
 * @param n：向量维度
 * @param dst：需要申请的向量
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,std::vector<ScalarType>& dst)
{
    dst.resize(n);
}



template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Fill：利用POD数据进行向量填充
 * @param src：POD数据源
 * @param n：向量维度
 * @param dst：目标向量
 */
CHIPSUM_FUNCTION_INLINE void Fill(
        const ScalarType* src,
        const SizeType n,
        std::vector<ScalarType>& dst
        )
{
    dst = std::vector<ScalarType>(src,src+n);
}



template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot：向量x和y的内积操作
 * @param x：向量x
 * @param y：向量y
 * @param n：向量维度
 * @param r：结果
 */
CHIPSUM_FUNCTION_INLINE void Dot(
        const std::vector<ScalarType>& x,

        const std::vector<ScalarType>& y,

        const SizeType& n,
        ScalarType& r
        )
{




    for(SizeType i=0;i<x.size();++i)
    {
        r += x[i]*y[i];
    }



}

template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Scal
 * @param R：输出，比例缩放后的向量
 * @param a：输入，缩放比例
 * @param X：输入，原向量
 */
CHIPSUM_FUNCTION_INLINE void Scal(std::vector<ScalarType>& R,const ScalarType a,const std::vector<ScalarType>& X){
    if(R.size() != X.size()) R.resize(X.size());
    for(size_t i=0;i<X.size();++i)
    {
        R[i] = a*X[i];
    }
}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Norm1：1范数
 * @param X：原向量
 * @return 范数结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm1(const std::vector<ScalarType>& X)
{
    ScalarType acc = 0.0;
    for(size_t i=0;i<X.size();++i){
        acc += X[i];
    }
    return acc;

}

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Norm2：2范数
 * @param X：原向量
 * @return 结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm2(const std::vector<ScalarType>& X)
{
    ScalarType acc = 0.0;
    for(size_t i=0;i<X.size();++i){
        acc += X[i]*X[i];
    }

    return std::sqrt(acc);

}



template<typename ScalarType,typename SizeType,typename ...Props>

/**
 * @brief Axpby： y(i) = a*x(i)+b*y(i)
 * @param a： 标量，如int和double
 * @param X： 向量X
 * @param Y： 向量Y
 */
CHIPSUM_FUNCTION_INLINE void Axpy(
        ScalarType a,
        const std::vector<ScalarType>& X,
        std::vector<ScalarType>& Y)
{
//    std::assert(X.size()==Y.size());
    for(size_t i=0;i<Y.size();++i)
    {
        Y[i] += a*X[i];
    }
}


template<typename ScalarType,typename SizeType,typename ...Props>

/**
 * @brief Axpby： y(i) = a*x(i)+b*y(i)
 * @param a： 标量，如int和double
 * @param X： 向量X
 * @param b： 标量，同a
 * @param Y： 向量Y
 */
CHIPSUM_FUNCTION_INLINE void Axpby(
        ScalarType a,
        std::vector<ScalarType>& X,
        ScalarType b,
        std::vector<ScalarType>& Y)
{
//    std::assert(X.size()==Y.size());
    for(size_t i=0;i<Y.size();++i)
    {
        Y[i] = a*X[i] + b*Y[i];
    }
}


template<typename ScalarType,typename SizeType,typename ...Props>

/**
 * @brief DeepCopy：深拷贝操作（该接口主要应对Kokkos等后端的拷贝机制）
 * @param dst：目标向量
 * @param src：原向量
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(
        std::vector<ScalarType>& dst,
        const std::vector<ScalarType>& src
        )
{

    for(size_t i=0;i<dst.size();++i)
    {
        dst[i] = src[i];
    }


}

template<typename ScalarType,typename SizeType,typename ...Props>

/**
 * @brief DeepCopy：浅拷贝操作（该接口主要应对Kokkos等后端的拷贝机制）
 * @param dst：目标向量
 * @param src：原向量
 */
CHIPSUM_FUNCTION_INLINE void ShallowCopy(
        std::vector<ScalarType>& dst,
        const std::vector<ScalarType>& src
        )
{
    dst = src;

}



template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Print：打印向量数据。
 * @param out：输出流，可以是std::cout，也可以是std::ofstream等
 * @param vec：向量
 */
CHIPSUM_FUNCTION_INLINE void Print(std::ostream& out,const std::vector<ScalarType>& vec)
{

    out<<" [";
    for(size_t i=0;i<vec.size()-1;++i)
    {
        out<<vec[i]<<", ";
    }

    out<<vec[vec.size()-1]<<"]"<<std::endl;


}



}// End namespace Vector
}// End namespace Impl

}// End namespace Numeric

}// End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
