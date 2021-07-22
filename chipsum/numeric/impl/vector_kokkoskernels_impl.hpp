#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__



#include <Kokkos_Core.hpp>
#include <KokkosBlas1_dot.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_fill.hpp>



#include "../numeric_traits.hpp"

static int vector_name = 0;
namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,ChipSum::Backend::KokkosKernels,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::KokkosKernels,Props...>
{
    using vector_type = typename Kokkos::View<ScalarType*>;


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
void CreateVector(const SizeType n,Kokkos::View<ScalarType*>& dst)
{
    dst = Kokkos::View<ScalarType*>("v"+std::to_string(vector_name),n);
}


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
        Kokkos::View<ScalarType*> dst
        )
{
    typename Kokkos::View<ScalarType*>::HostMirror h_dst("hst",dst.extent(0));

    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0,n),
                         KOKKOS_LAMBDA(const int i){

        h_dst(i) = src[i];

    }

    );

    Kokkos::deep_copy(dst,h_dst);

}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot 向量内积的串行实现
 * @param a 利用traits技术实现的向量类型
 * @param b
 * @param r
 */
void Dot(
        const Kokkos::View<ScalarType*>& a,

        const Kokkos::View<ScalarType*>& b,
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
