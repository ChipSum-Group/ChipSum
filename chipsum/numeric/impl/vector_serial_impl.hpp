#ifndef __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__

#include <vector>
#include <iostream>
using namespace std;

#include <Kokkos_Core.hpp>

#include "../numeric_traits.hpp"


namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,ChipSum::Backend::CPUSerialBackend,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUSerialBackend,Props...>
{
    using vector_type = typename std::vector<ScalarType>;


    using vector_type_reference = typename std::add_lvalue_reference<std::vector<ScalarType>>::type ;
    using const_vector_type_reference = typename std::add_const<vector_type_reference>::type;

    using scalar_type = typename std::vector<ScalarType>::value_type;
    using scalar_type_reference = typename std::add_lvalue_reference<typename std::vector<ScalarType>::value_type>::type;


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
void FillVector(
        ScalarType* src,
        SizeType n,
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
void Dot(
        const std::vector<ScalarType>& a,

        const std::vector<ScalarType>& b,

        const SizeType& n,
        ScalarType& r
        )
{

    cout<<"BBB"<<endl;

    r = 0.0;

    for(SizeType i=0;i<a.size();++i)
    {
        r += a[i]*b[i];
    }


}

}// End namespace Vector
}// End namespace Impl

}// End namespace Numeric

}// End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
