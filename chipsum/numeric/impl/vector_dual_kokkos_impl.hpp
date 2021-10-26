/*
 * @Description: 向量vector的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-09-06 09:02:26
 */

#ifndef __CHIPSUM_DUALVECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_DUALVECTOR_KOKKOSKERNELS_IMPL_HPP__

#include <fstream>
//#include <Kokkos_Core.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_scal.hpp>
#include <Kokkos_Vector.hpp>

#include <iostream>
using namespace std;

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"


namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Vector_Traits<ScalarType, SizeType, ChipSum::Backend::Kokkos,
                     Props...>
    : public Operator_Traits<ScalarType, SizeType,
                             ChipSum::Backend::Kokkos> {
  using vector_type = typename Kokkos::vector<ScalarType>;
  using size_type = std::size_t;

  using device_scalar_value_type = typename Kokkos::View<ScalarType>;
};

namespace Impl {

namespace Vector {


template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 函数子，用于响亮的scal操作（为满足Device端接口）
 * @author: Li Kunyun
 */
struct Vec_Scal_Functor {
  Vec_Scal_Functor(Kokkos::View<ScalarType> ai, Kokkos::vector<ScalarType> xi,
               Kokkos::vector<ScalarType> yi)
  {
    a = ai;
    x = xi;
    y = yi;
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    y.d_view(i) = a() * x.d_view(i);
  }

private:
  Kokkos::View<ScalarType> a;
  Kokkos::vector<ScalarType> x;
  Kokkos::vector<ScalarType> y;
};


//template <typename ScalarType, typename SizeType, typename... Props>
///**
// * @description: 函数子，用于响亮的scal操作（为满足Device端接口）
// * @author: Li Kunyun
// */
//struct Vec_Set_Functor_host {
//  Vec_Set_Functor(ScalarType* ai, Kokkos::vector<ScalarType> xi) {
//    a = ai;
//    x = xi;
//    y = yi;
//  }

//  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
//    y.d_view(i) = a() * x.d_view(i);
//  }

//private:
//  Kokkos::View<ScalarType> a;
//  Kokkos::vector<ScalarType> x;
//  Kokkos::vector<ScalarType> y;
//};


template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建未初始化的向量
 * @param {*} n 向量长度 
 * @param {*} x 向量（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,
                                    Kokkos::vector<ScalarType> &x) {

  x = Kokkos::vector<ScalarType>(static_cast<size_t>(n));

}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description: 创建初始化的向量
 * @param {*} src POD数据源
 * @param {*} n 向量长度
 * @param {*} x 向量（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Create(ScalarType *src, const std::size_t n,
                                    Kokkos::vector<ScalarType> &x) {



  x.set_overallocation(0.0);
  x.h_view = typename Kokkos::vector<ScalarType>::t_host(src,n);

  x.host_to_device();

}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 向量dot运算（by reference） = typename Kokkos::View<ScalarType*>::HostMirror(src,n);
 * @param {*} a 向量
 * @param {*} b 向量 = typename Kokkos::View<ScalarType*>::HostMirror(src,n);
 * @param {*} n 向量长度
 * @param {*} r 结果（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType Dot(const Kokkos::vector<ScalarType> &a,
                                 const Kokkos::vector<ScalarType> &b,
                                 const SizeType n) {

  return KokkosBlas::dot(a.d_view, b.d_view);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 向量dot运算（by reference）
 * @param {*} a 向量
 * @param {*} b 向量
 * @param {*} n 向量长度
 * @param {*} r 结果（out）
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void
Dot(const Kokkos::vector<ScalarType> &x, const Kokkos::vector<ScalarType> &y,
    const SizeType n, Kokkos::View<ScalarType> &r) {

  KokkosBlas::dot(r, x.d_view, y.d_view);
}



template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: R=a*X
 * @param {*} R Scal结果（out）
 * @param {*} a X的系数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Scal(Kokkos::vector<ScalarType> &R,
                                  const Kokkos::View<ScalarType> &a,
                                  const Kokkos::vector<ScalarType> &X) {
  Kokkos::parallel_for(R.extent(0),
                       Vec_Scal_Functor<ScalarType, SizeType>(a, X, R));

  R.on_device();
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description: R=a*X
 * @param {*} R Scal结果（out）
 * @param {*} a X的系数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Scal(Kokkos::vector<ScalarType> &R,
                                  const ScalarType &a,
                                  const Kokkos::vector<ScalarType> &X) {
  KokkosBlas::scal(R.d_view, a, X.d_view);
  R.on_device();
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: X的1范数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm1(const Kokkos::vector<ScalarType> &X) {
  return KokkosBlas::nrm1(X.d_view);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: X的2范数
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm2(const Kokkos::vector<ScalarType> &X) {
  return KokkosBlas::nrm2(X.d_view);
}


template <typename ScalarType, typename SizeType, typename Arg,typename... Props>
/**
 * @description: X的1范数
 * @param {*} r 结果
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Norm1(Arg& r,const Kokkos::vector<ScalarType> &X) {
  KokkosBlas::nrm1(r,X.d_view);
}

template <typename ScalarType, typename SizeType, typename Arg,typename... Props>
/**
 * @description: X的2范数
 * @param {*} r 结果
 * @param {*} X 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Norm2(Arg& r,const Kokkos::vector<ScalarType> &X) {
  KokkosBlas::nrm2(r,X.d_view);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: Y=Y+a*X
 * @param {*} a 系数
 * @param {*} X 向量
 * @param {*} Y 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Axpy(ScalarType a,
                                  const Kokkos::vector<ScalarType> &X,
                                  const Kokkos::vector<ScalarType> &Y) {
  KokkosBlas::axpy(a, X.d_view, Y.d_view);
  Y.on_device();
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: Y=b*Y+a*X
 * @param {*} a 系数
 * @param {*} X 向量
 * @param {*} b 系数
 * @param {*} Y 向量
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void
Axpby(ScalarType a, const Kokkos::vector<ScalarType> &X, ScalarType b,
      const Kokkos::vector<ScalarType> &Y) {
  KokkosBlas::axpby(a, X, b, Y);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 深拷贝
 * @param {*} dst 目标数据
 * @param {*} src 原数据
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(const Kokkos::vector<ScalarType> &dst,
                                      const Kokkos::vector<ScalarType> &src) {

  Kokkos::deep_copy(dst, src);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 浅拷贝
 * @param {*} dst 目标数据
 * @param {*} src 原数据
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void
ShallowCopy(const Kokkos::vector<ScalarType> &dst,
            const Kokkos::vector<ScalarType> &src) {
  dst = src;
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 获取向量元素（GPU端该函数会报错）
 * @param {*} index 索引
 * @param {*} vec 向量
 * @return {*} 元素值
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE ScalarType &GetItem(const std::size_t index,
                                            Kokkos::vector<ScalarType> &vec) {

  return vec(index);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印向量，一般用于调试
 * @param {*} vec 向量
 * @param {*} out 输出流
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Print(const Kokkos::vector<ScalarType> &vec,
  std::ostream &out
                                   ) {
  typename Kokkos::vector<ScalarType>::HostMirror h_vec("h_vector",
                                                        vec.extent(0));
  Kokkos::deep_copy(h_vec, vec);

  out << vec.label() << ": [";
  for (std::size_t i = 0; i < h_vec.extent(0) - 1; ++i) {
    out << h_vec(i) << ", ";
  }

  out << h_vec(h_vec.extent(0) - 1) << "]" << std::endl;
}



} // End namespace Vector

} // End namespace Impl

} // End namespace Numeric

} // End namespace ChipSum

#endif // __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
