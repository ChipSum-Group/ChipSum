/* * * * * * * * * * * * * * * * * * * * *
*   File:     vector_kokkoskernels_impl.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__


#include <fstream>
//#include <Kokkos_Core.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <Kokkos_Vector.hpp>



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

    using device_scalar_value_type = typename Kokkos::View<ScalarType>;

};


namespace Impl {

namespace Vector {

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Create：创建向量数据，主要是申请内存空间
 * @param n：向量维度
 * @param dst：需要申请的向量
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,Kokkos::View<ScalarType*>& dst)
{

    dst = Kokkos::View<ScalarType*>("vector_"+std::to_string(vector_name++),static_cast<size_t>(n) );
}



template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Fill：利用POD数据进行向量填充
 * @param src：POD数据源
 * @param n：向量维度
 * @param dst：目标向量
 */
CHIPSUM_FUNCTION_INLINE void Create(
        ScalarType* src,
        const std::size_t n,
        Kokkos::View<ScalarType*>& dst
        )
{
    typename Kokkos::View<ScalarType*>::HostMirror h_dst(src,n);

    if(dst.extent(0) != n) {
        dst = Kokkos::View<ScalarType*>("vector_"+std::to_string(vector_name++),n);
    }
    Kokkos::deep_copy(dst,h_dst);




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
        const Kokkos::View<ScalarType*>& a,
        const Kokkos::View<ScalarType*>& b,
        const SizeType n,
        ScalarType& r
        )
{

    r = KokkosBlas::dot(a,b);

}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Dot：向量x和y的内积操作（未测试）
 * @param x：向量x
 * @param y：向量y
 * @param n：向量维度
 * @param r：结果（Kokkos标量数据）
 */
CHIPSUM_FUNCTION_INLINE void Dot(
        const Kokkos::View<ScalarType*>& a,
        const Kokkos::View<ScalarType*>& b,
        const SizeType n,
        Kokkos::View<ScalarType>& r
        )
{

    KokkosBlas::dot(r,a,b);
}



template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Scal
 * @param R：输出，比例缩放后的向量
 * @param a：输入，缩放比例
 * @param X：输入，原向量
 */
CHIPSUM_FUNCTION_INLINE void Scal(Kokkos::View<ScalarType*>& R,
                                  const ScalarType a,
                                  const Kokkos::View<ScalarType*>& X)
{
    KokkosBlas::scal(R,a,X);
}




template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Norm1：1范数
 * @param X：原向量
 * @return 范数结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm1(
        const Kokkos::View<ScalarType*>& X)
{
    return KokkosBlas::nrm1(X);
}




template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Norm2：2范数
 * @param X：原向量
 * @return 结果
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm2(
        const Kokkos::View<ScalarType*>& X)
{
    return KokkosBlas::nrm2(X);
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
        const Kokkos::View<ScalarType*>& X,
        const Kokkos::View<ScalarType*>& Y)
{
    KokkosBlas::axpy(a,X,Y);
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
        const Kokkos::View<ScalarType*>& X,
        ScalarType b,
        const Kokkos::View<ScalarType*>& Y)
{
    KokkosBlas::axpby(a,X,b,Y);
}


template<typename ScalarType,typename SizeType,typename ...Props>

/**
 * @brief DeepCopy：深拷贝操作（该接口主要应对Kokkos等后端的拷贝机制）
 * @param dst：目标向量
 * @param src：原向量
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(
        const Kokkos::View<ScalarType*>& dst,
        const Kokkos::View<ScalarType*>& src
        )
{

    Kokkos::deep_copy(dst,src);

}

template<typename ScalarType,typename SizeType,typename ...Props>

/**
 * @brief DeepCopy：浅拷贝操作（该接口主要应对Kokkos等后端的拷贝机制）
 * @param dst：目标向量
 * @param src：原向量
 */
CHIPSUM_FUNCTION_INLINE void ShallowCopy(
        const Kokkos::View<ScalarType*>& dst,
        const Kokkos::View<ScalarType*>& src
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
CHIPSUM_FUNCTION_INLINE void Print(
        std::ostream& out,
        const Kokkos::View<ScalarType*>& vec)
{
    typename Kokkos::View<ScalarType*>::HostMirror h_vec("h_vector",vec.extent(0));
    Kokkos::deep_copy(h_vec,vec);

    out<<vec.label()<<": [";
    for(size_t i=0;i<h_vec.extent(0)-1;++i)
    {
        out<<h_vec(i)<<", ";
    }

    out<<h_vec(h_vec.extent(0)-1)<<"]"<<std::endl;


}

template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief GetItem：获取下标对应的值
 * @param index：下标
 * @param vec：原向量
 * @return
 */
CHIPSUM_FUNCTION_INLINE ScalarType& GetItem(const std::size_t index,Kokkos::View<ScalarType*>& vec)
{

    return vec(index);

}

}// End namespace Vector

}// End namespace Impl

}// End namespace Numeric

}// End namespace ChipSum

#endif // __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
