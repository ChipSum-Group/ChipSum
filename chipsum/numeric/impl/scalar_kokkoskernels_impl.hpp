/*
 * @Author: your name
 * @Date: 2021-08-09 12:34:28
 * @LastEditTime: 2021-08-10 11:02:00
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: scalar_kokkoskernels_impl.hpp
 */

#ifndef __CHIPSUM_SCALAR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_SCALAR_KOKKOSKERNELS_IMPL_HPP__

#include <KokkosBlas1_axpby_spec.hpp>



#include "../numeric_traits.hpp"
#include "../../backend/backend.hpp"
#include "../../chipsum_macro.h"

static int scalar_name = 0;

namespace ChipSum
{
    namespace Numeric
    {

        template <typename ScalarType, typename SizeType, typename... Props>
        struct Scalar_Traits<ScalarType, SizeType, ChipSum::Backend::KokkosKernels, Props...> : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::KokkosKernels, Props...>
        {
            using scalar_type = Kokkos::View<ScalarType>;

            using device_scalar_value_type = typename scalar_type::HostMirror;
        };

        namespace Impl
        {
            namespace Scalar
            {

                template <typename ScalarType, typename SizeType, typename... Props>
                CHIPSUM_FUNCTION_INLINE void Create(Kokkos::View<ScalarType> &r)
                {
                    r = Kokkos::View<ScalarType>("scalar_" + std::to_string(scalar_name++));
                }

                template <typename ScalarType, typename SizeType, typename... Props>
                CHIPSUM_FUNCTION_INLINE void DeepCopy(const ScalarType s, Kokkos::View<ScalarType> &r)
                {

                    typename Kokkos::View<ScalarType>::HostMirror h_r = Kokkos::create_mirror_view(r);
                    h_r() = s;
                    Kokkos::deep_copy(r, h_r);
                }

                template <typename ScalarType, typename SizeType, typename... Props>
                CHIPSUM_FUNCTION_INLINE void Create(const ScalarType s, Kokkos::View<ScalarType> &r)
                {
                    r = Kokkos::View<ScalarType>("scalar_" + std::to_string(scalar_name++));
                    DeepCopy<ScalarType, SizeType>(s, r);
                }

                template <typename ScalarType, typename SizeType, typename... Props>
                CHIPSUM_FUNCTION_INLINE const ScalarType GetItem(const Kokkos::View<ScalarType> &s)
                {
                    typename Kokkos::View<ScalarType>::HostMirror h_s = Kokkos::create_mirror_view(s);
                    Kokkos::deep_copy(h_s, s);
                    return h_s();
                }

                template <typename ScalarType, typename SizeType, typename... Props>
                CHIPSUM_FUNCTION_INLINE void GetItem(const Kokkos::View<ScalarType> &s, ScalarType &r)
                {
                    typename Kokkos::View<ScalarType>::HostMirror h_s = Kokkos::create_mirror_view(s);
                    Kokkos::deep_copy(h_s, s);
                    r = h_s();
                }

                template <typename ScalarType, typename SizeType, typename... Props>
                CHIPSUM_FUNCTION_INLINE void Mult(const Kokkos::View<ScalarType> &s,
                                                  const Kokkos::View<ScalarType *> &v,
                                                  Kokkos::View<ScalarType *> &r)
                {
                    
                    KokkosBlas::scal(r,s(),v);
                }

            }
        }
    }
}

#endif //
