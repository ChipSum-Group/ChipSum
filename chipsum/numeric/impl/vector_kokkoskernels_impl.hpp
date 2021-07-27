#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__



#include <Kokkos_Core.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_axpby.hpp>




#include "../numeric_traits.hpp"
#include "../../chipsum_macro.h"

static int vector_name = 0;
namespace ChipSum{
namespace Numeric {


template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,ChipSum::Backend::KokkosKernels,Props...>
        :
        public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::KokkosKernels>
{
    using vector_type = typename Kokkos::View<ScalarType*>;
    using size_type = std::size_t;

};


namespace Impl {

namespace Vector {

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief 创建向量
 * @param n：向量维度
 * @param dst：向量
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,Kokkos::View<ScalarType*>& dst)
{
    dst = Kokkos::View<ScalarType*>("vector_"+std::to_string(vector_name++),n);
}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief POD类型数组填充向量
 * @param POD数据
 * @param n
 * @param dst
 */
CHIPSUM_FUNCTION_INLINE void Fill(
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
CHIPSUM_FUNCTION_INLINE void Dot(
        const Kokkos::View<ScalarType*>& a,

        const Kokkos::View<ScalarType*>& b,
        const SizeType n,
        ScalarType& r
        )
{
    r = KokkosBlas::dot(a,b);
}

template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Scal
 * @param R：输出，比例缩放后的向量
 * @param a：输入，缩放比例
 * @param X：输入，原向量
 */
CHIPSUM_FUNCTION_INLINE void Scal(Kokkos::View<ScalarType*>& R,const ScalarType a,const Kokkos::View<ScalarType*>& X){
    KokkosBlas::scal(R,a,X);
}

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Norm1：1范数
 * @param X：原向量
 * @return 范数结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm1(const Kokkos::View<ScalarType*>& X)
{
    return KokkosBlas::nrm1(X);
}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Norm2：2范数
 * @param X：原向量
 * @return 范数结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm2(const Kokkos::View<ScalarType*>& X)
{
    return KokkosBlas::nrm2(X);
}

}// End namespace Vector

}// End namespace Impl

}// End namespace Numeric

}// End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
