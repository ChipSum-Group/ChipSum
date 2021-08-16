/*
 * @Author: your name
 * @Date: 2021-08-09 12:34:28
 * @LastEditTime: 2021-08-16 10:24:47
 * @LastEditors: Li Kunyun
 * @Description: In User Settings Edit
 * @FilePath: scalar_kokkoskernels_impl.hpp
 */

#ifndef __CHIPSUM_SCALAR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_SCALAR_KOKKOSKERNELS_IMPL_HPP__

#include <KokkosBlas1_axpby_spec.hpp>
#include <fstream>

#include "../../backend/backend.hpp"
#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

static int scalar_name = 0;

namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Scalar_Traits<ScalarType, SizeType, ChipSum::Backend::KokkosKernels,
                     Props...>
    : public Operator_Traits<ScalarType, SizeType,
                             ChipSum::Backend::KokkosKernels, Props...> {
  using scalar_type = Kokkos::View<ScalarType>;

  using device_scalar_value_type = typename scalar_type::HostMirror;
};

namespace Impl {
namespace Scalar {

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建标量
 * @param {View<ScalarType>} &r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Create(Kokkos::View<ScalarType> &r) {
  r = Kokkos::View<ScalarType>("scalar_" + std::to_string(scalar_name++));
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 标量深拷贝
 * @param {*} s 数据源
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void DeepCopy(const ScalarType s,
                                      Kokkos::View<ScalarType> &r) {

  typename Kokkos::View<ScalarType>::HostMirror h_r =
      Kokkos::create_mirror_view(r);
  h_r() = s;
  Kokkos::deep_copy(r, h_r);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建标量
 * @param {*} s 数据源
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Create(const ScalarType s,
                                    Kokkos::View<ScalarType> &r) {
  DeepCopy<ScalarType, SizeType>(s, r);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 获取标量（by return）
 * @param {*} s 标量
 * @return {*} 标量（out）
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE const ScalarType
GetItem(const Kokkos::View<ScalarType> &s) {
  typename Kokkos::View<ScalarType>::HostMirror h_s =
      Kokkos::create_mirror_view(s);
  Kokkos::deep_copy(h_s, s);
  return h_s();
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 获取标量（by reference）
 * @param {*} s 标量
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void GetItem(const Kokkos::View<ScalarType> &s,
                                     ScalarType &r) {
  typename Kokkos::View<ScalarType>::HostMirror h_s =
      Kokkos::create_mirror_view(s);
  Kokkos::deep_copy(h_s, s);
  r = h_s();
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印标量，一般用于调试
 * @param {*} s 标量
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Print(const Kokkos::View<ScalarType> &s,
                                   std::ostream &out) {
  typename Kokkos::View<ScalarType>::HostMirror h_s =
      Kokkos::create_mirror_view(s);
  Kokkos::deep_copy(h_s, s);
  out << s.label() << ": ";
  out << h_s() << endl;
}





} // namespace Scalar
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif //
