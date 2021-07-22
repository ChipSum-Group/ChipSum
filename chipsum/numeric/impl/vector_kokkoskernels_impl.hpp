#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__


#include <list>
#include "../numeric_traits.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_dot.hpp>


namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,ChipSum::Backend::KokkosKernels,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::KokkosKernels,Props...>
{
    using vector_type = typename Kokkos::View<ScalarType*,Kokkos::DefaultExecutionSpace>;


    using vector_type_reference = typename std::add_lvalue_reference<vector_type>::type ;
    using const_vector_type_reference = typename std::add_const<vector_type_reference>::type;

    using scalar_type = ScalarType;
    using scalar_type_reference = typename std::add_lvalue_reference<ScalarType>::type;


    using size_type = SizeType;

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
        const SizeType n,
        std::list<ScalarType>& dst
        )
{
    dst = std::list<ScalarType>(src,src+n);
}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot 向量内积的串行实现
 * @param a 利用traits技术实现的向量类型
 * @param b
 * @param r
 */
void Dot(
        const Kokkos::View<ScalarType*,Kokkos::DefaultExecutionSpace>& a,

        const Kokkos::View<ScalarType*,Kokkos::DefaultExecutionSpace>& b,
        const SizeType n,
        ScalarType& r
        )
{




    r = KokkosBlas::dot(a,b);


}


}// End namespace Vector

}// End namespace Impl

}// End namespace Numeric

}// End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
