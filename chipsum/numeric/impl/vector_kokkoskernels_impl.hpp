/*
 * @Description: 向量vector的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 16:17:59
 */

#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__

#include <fstream>
//#include <Kokkos_Core.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <KokkosBlas1_scal.hpp>
#include <Kokkos_Vector.hpp>


#include <iostream>
using namespace std;

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"


/// 一个不幸的消息是vector的实现需要大改。
/// 由于直接用view<type*>作为vector的数据底层，导致了数据无法灵活在
/// device和host之间流转。所以这部分的底层数据类型需要替换为dualview
/// 这个改动将会带来整个kokkos实现的大改，但好在接口本就不多。

/// 这部分工作不会太复杂。

static int vector_name = 0;
namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Vector_Traits<ScalarType, SizeType, ChipSum::Backend::KokkosKernels,
                     Props...>
    : public Operator_Traits<ScalarType, SizeType,
                             ChipSum::Backend::KokkosKernels> {
  using vector_type = typename Kokkos::View<ScalarType *>;
  using size_type = std::size_t;

  using device_scalar_value_type = typename Kokkos::View<ScalarType>;
};

namespace Impl {

namespace Vector {

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void create(const SizeType n,
                                    Kokkos::View<ScalarType *> &x) {

  x = Kokkos::View<ScalarType *>("vector_" + std::to_string(vector_name++),
                                   static_cast<size_t>(n));
}

template <typename ScalarType, typename SizeType, typename... Props>


CHIPSUM_FUNCTION_INLINE void create(ScalarType *src, const std::size_t n,
                                    Kokkos::View<ScalarType *> &x) {
  typename Kokkos::View<ScalarType *>::HostMirror h_x(src, n);

  if(x.extent(0)==0){
    x = Kokkos::View<ScalarType*>("vector_"+std::to_string(vector_name++),n);
  }
  Kokkos::deep_copy(x, h_x);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType dot(const Kokkos::View<ScalarType *> &a,
                                 const Kokkos::View<ScalarType *> &b,
                                 const SizeType n) {

  return KokkosBlas::dot(a, b);
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void
dot(const Kokkos::View<ScalarType *> &x, const Kokkos::View<ScalarType *> &y,
    const SizeType n, Kokkos::View<ScalarType> &r) {

  KokkosBlas::dot(r, x, y);
}

template <typename ScalarType, typename SizeType, typename... Props>
/// x = a*x 的函数子
struct scal_functor {
  scal_functor(Kokkos::View<ScalarType> ai, Kokkos::View<ScalarType *> xi,
               Kokkos::View<ScalarType *> yi) {
    a = ai;
    x = xi;
    y = yi;
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    y(i) = a() * x(i);
  }

private:
  Kokkos::View<ScalarType> a;
  Kokkos::View<ScalarType *> x;
  Kokkos::View<ScalarType *> y;
};

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void scal(Kokkos::View<ScalarType *> &R,
                                  const Kokkos::View<ScalarType> &a,
                                  const Kokkos::View<ScalarType *> &X) {
  Kokkos::parallel_for(R.extent(0),
                       scal_functor<ScalarType, SizeType>(a, X, R));
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void scal(Kokkos::View<ScalarType *> &R,
                                  const ScalarType &a,
                                  const Kokkos::View<ScalarType *> &X) {
  KokkosBlas::scal(R, a, X);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType norm1(const Kokkos::View<ScalarType *> &X) {
  return KokkosBlas::nrm1(X);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType norm2(const Kokkos::View<ScalarType *> &X) {
  return KokkosBlas::nrm2(X);
}




template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType norminf(const Kokkos::View<ScalarType *> &X) {
  return KokkosBlas::nrminf(X);
}





template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void axpy(ScalarType a,
                                  const Kokkos::View<ScalarType *> &X,
                                  const Kokkos::View<ScalarType *> &Y) {
  KokkosBlas::axpy(a, X, Y);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
axpby(ScalarType a, const Kokkos::View<ScalarType *> &X, ScalarType b,
      const Kokkos::View<ScalarType *> &Y) {
  KokkosBlas::axpby(a, X, b, Y);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void deep_copy(const Kokkos::View<ScalarType *> &dst,
                                      const Kokkos::View<ScalarType *> &src) {

  Kokkos::deep_copy(dst, src);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
shallow_copy(const Kokkos::View<ScalarType *> &dst,
            const Kokkos::View<ScalarType *> &src) {
  dst = src;
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE ScalarType &get_item(const std::size_t index,
                                            const Kokkos::View<ScalarType *> &vec) {

  return vec(index);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void print(const Kokkos::View<ScalarType *> &vec,
  std::ostream &out
                                   ) {
  typename Kokkos::View<ScalarType *>::HostMirror h_vec("h_vector",
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
