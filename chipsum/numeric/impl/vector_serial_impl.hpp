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
using namespace std;



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

};




namespace Impl {

namespace Vector {

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot 向量内积的串行实现
 * @param a 利用traits技术实现的向量类型
 * @param b
 * @param r
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,std::vector<ScalarType>& dst)
{
    dst.resize(n);
}



template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot 向量内积的串行实现
 * @param a 利用traits技术实现的向量类型
 * @param b
 * @param r
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
 * @brief Dot 向量内积的串行实现
 * @param a 利用traits技术实现的向量类型
 * @param b
 * @param r
 */
CHIPSUM_FUNCTION_INLINE void Dot(
        const std::vector<ScalarType>& a,

        const std::vector<ScalarType>& b,

        const SizeType& n,
        ScalarType& r
        )
{


    r = 0.0;

    for(SizeType i=0;i<a.size();++i)
    {
        r += a[i]*b[i];
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
 * @brief Norm1：1范数
 * @param X：原向量
 * @return 范数结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm2(const std::vector<ScalarType>& X)
{
    ScalarType acc = 0.0;
    for(size_t i=0;i<X.size();++i){
        acc += X[i]*X[i];
    }

    return std::sqrt(acc);

}

}// End namespace Vector
}// End namespace Impl

}// End namespace Numeric

}// End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
