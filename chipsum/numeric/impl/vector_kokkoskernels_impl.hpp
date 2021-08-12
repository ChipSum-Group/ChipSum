/*
 * @Description: 向量vector的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 14:18:23
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
#include <KokkosBlas1_scal.hpp>
#include <Kokkos_Vector.hpp>

#include <iostream>
using namespace std;

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

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
/**
 * @description:
 * @param {const SizeType} n
 * @param {View<ScalarType *>} &dst
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Create(const SizeType n,
                                    Kokkos::View<ScalarType *> &dst) {

  dst = Kokkos::View<ScalarType *>("vector_" + std::to_string(vector_name++),
                                   static_cast<size_t>(n));
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Create(ScalarType *src, const std::size_t n,
                                    Kokkos::View<ScalarType *> &dst) {
  typename Kokkos::View<ScalarType *>::HostMirror h_dst(src, n);

  if (dst.extent(0) != n) {
    dst = Kokkos::View<ScalarType *>("vector_" + std::to_string(vector_name++),
                                     n);
  }
  Kokkos::deep_copy(dst, h_dst);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Dot(const Kokkos::View<ScalarType *> &a,
                                 const Kokkos::View<ScalarType *> &b,
                                 const SizeType n, ScalarType &r) {

  r = KokkosBlas::dot(a, b);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void
Dot(const Kokkos::View<ScalarType *> &a, const Kokkos::View<ScalarType *> &b,
    const SizeType n, Kokkos::View<ScalarType> &r) {

  KokkosBlas::dot(r, a, b);
}

template <typename ScalarType, typename SizeType, typename... Props>
struct Scal_Functor {
  Scal_Functor(Kokkos::View<ScalarType> ai, Kokkos::View<ScalarType *> xi,
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

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Scal(Kokkos::View<ScalarType *> &R,
                                  const Kokkos::View<ScalarType> &a,
                                  const Kokkos::View<ScalarType *> &X) {
  assert(X.extent(0) == R.extent(0));
  Kokkos::parallel_for(R.extent(0),
                       Scal_Functor<ScalarType, SizeType>(a, X, R));
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Scal(Kokkos::View<ScalarType *> &R,
                                  const ScalarType &a,
                                  const Kokkos::View<ScalarType *> &X) {
  KokkosBlas::scal(R, a, X);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType Norm1(const Kokkos::View<ScalarType *> &X) {
  return KokkosBlas::nrm1(X);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */

CHIPSUM_FUNCTION_INLINE ScalarType Norm2(const Kokkos::View<ScalarType *> &X) {
  return KokkosBlas::nrm2(X);
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Axpy(ScalarType a,
                                  const Kokkos::View<ScalarType *> &X,
                                  const Kokkos::View<ScalarType *> &Y) {
  KokkosBlas::axpy(a, X, Y);
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void
Axpby(ScalarType a, const Kokkos::View<ScalarType *> &X, ScalarType b,
      const Kokkos::View<ScalarType *> &Y) {
  KokkosBlas::axpby(a, X, b, Y);
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(const Kokkos::View<ScalarType *> &dst,
                                      const Kokkos::View<ScalarType *> &src) {

  Kokkos::deep_copy(dst, src);
}

template <typename ScalarType, typename SizeType, typename... Props>

/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void
ShallowCopy(const Kokkos::View<ScalarType *> &dst,
            const Kokkos::View<ScalarType *> &src) {
  dst = src;
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description:
 * @param {*}
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out,
                                   const Kokkos::View<ScalarType *> &vec) {
  typename Kokkos::View<ScalarType *>::HostMirror h_vec("h_vector",
                                                        vec.extent(0));
  Kokkos::deep_copy(h_vec, vec);

  out << vec.label() << ": [";
  for (size_t i = 0; i < h_vec.extent(0) - 1; ++i) {
    out << h_vec(i) << ", ";
  }

  out << h_vec(h_vec.extent(0) - 1) << "]" << std::endl;
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description:
 * @param {const} std
 * @param {View<ScalarType *>} &vec
 * @return {*}
 */
CHIPSUM_FUNCTION_INLINE ScalarType &GetItem(const std::size_t index,
                                            Kokkos::View<ScalarType *> &vec) {

  return vec(index);
}

} // End namespace Vector

} // End namespace Impl

} // End namespace Numeric

} // End namespace ChipSum

#endif // __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
