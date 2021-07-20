#ifndef __VECTOR_SERIAL_IMPL_HPP__
#define __VECTOR_SERIAL_IMPL_HPP__

#include <vector>
#include "../numeric_traits.hpp"

#include <iostream>
using namespace std;

namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,ChipSum::Backend::CPUSerialBackend,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUSerialBackend,Props...>
{
    using vector_type = typename std::vector<ScalarType>;
//    using vector_type = typename std::vector<ScalarType>;


    using vector_type_reference = typename std::add_lvalue_reference<std::vector<ScalarType>>::type ;
    using const_vector_type_reference = typename std::add_const<vector_type_reference>::type;
    using scalar_type = typename std::vector<ScalarType>::value_type;
    using scalar_type_reference = typename std::add_lvalue_reference<typename std::vector<ScalarType>::value_type>::type;
    using size_type = typename std::vector<ScalarType>::size_type;

};










namespace Impl {


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot 向量内积的串行实现
 * @param a 利用traits技术实现的向量类型
 * @param b
 * @param r
 */
void VectorCopy(
        ScalarType* src,
        typename
        Vector_Traits<ScalarType,SizeType,
            ChipSum::Backend::CPUSerialBackend,Props...>::size_type n,
        typename
        Vector_Traits<ScalarType,SizeType,
            ChipSum::Backend::CPUSerialBackend,Props...>::vector_type_reference dst
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
void Dot(
        typename
        Vector_Traits<ScalarType,SizeType,
            ChipSum::Backend::CPUSerialBackend,Props...>::const_vector_type_reference a,

        typename
        Vector_Traits<ScalarType,SizeType,
                ChipSum::Backend::CPUSerialBackend,Props...>::const_vector_type_reference b,
        typename
        Vector_Traits<ScalarType,SizeType,
                ChipSum::Backend::CPUSerialBackend,Props...>::scalar_type_reference r
        )
{


    using traits = Vector_Traits<ScalarType,SizeType,
    ChipSum::Backend::CPUSerialBackend,Props...>;


    static_assert (std::is_same<decltype (a),decltype (b)>::value,"Parameter a and b shall be same type." );

    r = 0.0;
    typename traits::size_type n = a.size();

    for(typename traits::size_type i=0;i<n;++i)
    {
        r += a[i]*b[i];
    }


}


}

}

}

#endif
