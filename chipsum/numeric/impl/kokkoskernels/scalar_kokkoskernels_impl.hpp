///
/// \file     scalar_kokkoskernels_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-25
/// \brief    %stuff%
///

#ifndef __CHIPSUM_SCALAR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_SCALAR_KOKKOSKERNELS_IMPL_HPP__

#include <KokkosBlas1_axpby_spec.hpp>
#include <fstream>

#include "../../../backend/backend.hpp"
#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"

static int scalar_name = 0;

namespace ChipSum {
namespace Numeric {

template <typename ValueType>
struct Scalar_Traits<ValueType,ChipSum::Backend::KokkosKernels>
        : public Operator_Traits<ValueType> {
    using scalar_type = typename Kokkos::View<ValueType>;

    using host_type = typename scalar_type::HostMirror;

    using value_type = typename scalar_type::value_type;
};



namespace Impl {
namespace Scalar {

template<typename ValueType> using traits = Scalar_Traits<ValueType,ChipSum::Backend::KokkosKernels>;

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void create(Kokkos::View<ValueType>& r) {
    r = Kokkos::View<ValueType>("scalar_" + std::to_string(scalar_name++));
}


template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void deep_copy(
        Kokkos::View<ValueType>& r,
        const ValueType& s
                                       ) {

    typename Kokkos::View<ValueType>::HostMirror h_r
            = Kokkos::create_mirror_view(r);
    h_r() = s;
    Kokkos::deep_copy(r, h_r);
}

template <typename ValueType>
/**
 * @description: 创建标量
 * @param {*} s 数据源
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void create(
        const Kokkos::View<ValueType>& r,
        const ValueType& s
        ) {
    deep_copy(s, r);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE ValueType&
get_item(Kokkos::View<ValueType>& s) {
    typename Kokkos::View<ValueType>::HostMirror h_s
            = Kokkos::create_mirror_view(s);
    Kokkos::deep_copy(h_s, s);
    return h_s();
}

template <typename ValueType>
/**
 * @description: 获取标量（by reference）
 * @param {*} s 标量
 * @param {*} r 标量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
get_item(const Kokkos::View<ValueType>& s,
         ValueType& r)
{
    typename Kokkos::View<ValueType>::HostMirror h_s
            = Kokkos::create_mirror_view(s);
    Kokkos::deep_copy(h_s, s);
    r = h_s();
}

template <typename ValueType, typename OStreamT>
CHIPSUM_FUNCTION_INLINE
void print(const Kokkos::View<ValueType>& s,
           OStreamT& out) {
    typename Kokkos::View<ValueType>::HostMirror h_s =
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
