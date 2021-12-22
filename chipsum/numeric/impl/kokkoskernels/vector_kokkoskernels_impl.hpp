///
/// \file     vector_kokkoskernels_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-26
/// \brief    %stuff%
///

#ifndef __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__



#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)

#include <fstream>

#include <KokkosBlas1_fill.hpp>




#include <iostream>
using namespace std;

#include "../impl_abstract.hpp"
#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"

#include "vector_kokkoskernels_dot_impl.hpp"
#include "vector_kokkoskernels_scal_impl.hpp"
#include "vector_kokkoskernels_nrm1_impl.hpp"
#include "vector_kokkoskernels_nrm2_impl.hpp"
#include "vector_kokkoskernels_nrminf_impl.hpp"
#include "vector_kokkoskernels_axpby_impl.hpp"

/// 一个不幸的消息是vector的实现需要大改。
/// 由于直接用view<type*>作为vector的数据底层，导致了数据无法灵活在
/// device和host之间流转。所以这部分的底层数据类型需要替换为dualview
/// 这个改动将会带来整个kokkos实现的大改，但好在接口本就不多。

/// 这部分工作不会太复杂。

static int vector_name = 0;
namespace ChipSum {
namespace Numeric {

template <typename ValueType,typename ...Props>
struct Vector_Traits<ValueType, ChipSum::Backend::KokkosKernels,Props...>
        : public Operator_Traits<ValueType> {
    using vector_type = Kokkos::View<ValueType*>;
    using size_type = typename vector_type::size_type;
    using scalar_type = Kokkos::View<ValueType>;
    using value_type = typename vector_type::value_type;

    using vector_type_ref  =
    typename ::std::add_lvalue_reference<vector_type>::type;
    using const_vector_type_ref =
    typename ::std::add_const<vector_type_ref>::type;


    using size_type_ref =
    typename ::std::add_lvalue_reference<size_type>::type;
    using const_size_type_ref =
    typename ::std::add_const<size_type_ref>::type;



    using scalar_type_ref =
    typename ::std::add_lvalue_reference<scalar_type>::type;
    using const_scalar_type_ref =
    typename ::std::add_const<scalar_type_ref>::type;

    using value_type_ref =
    typename ::std::add_lvalue_reference<value_type>::type;
    using const_value_type_ref =
    typename ::std::add_const<value_type_ref>::type;

};



namespace Impl {
namespace Vector {


template <typename ValueType> using traits = Vector_Traits<ValueType,ChipSum::Backend::KokkosKernels>;

template<typename ValueType,typename ST>

CHIPSUM_FUNCTION_INLINE void create(
        Kokkos::View<ValueType*>& x,
        const ST& n
        )
{

    x = typename traits<ValueType>::vector_type("vector_" + ::std::to_string(vector_name++),
                                                n);
}

template <typename ValueType,typename ST>
CHIPSUM_FUNCTION_INLINE void create(
        Kokkos::View<ValueType*>& x,
        const ST& n,
        ValueType* src
        ) {
    typename traits<ValueType>::vector_type::HostMirror h_x(src, n);

    if(x.extent(0)!=h_x.extent(0)){
        x = typename traits<ValueType>::vector_type("vector_"+ ::std::to_string(vector_name++),n);
    }
    Kokkos::deep_copy(x, h_x);
}



template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
deep_copy(
        const Kokkos::View<ValueType*>& dst,
        const Kokkos::View<ValueType*>& src)
{

    Kokkos::deep_copy(dst, src);
}

template <typename ValueType,typename ST>
CHIPSUM_FUNCTION_INLINE
typename ::std::add_lvalue_reference<
typename traits<ValueType>::vector_type::value_type
>::type
get_item(
        const Kokkos::View<ValueType*>& vec,
        const ST& index
        )
{
    return vec(index);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE
void print(const Kokkos::View<ValueType*>& vec,
           ::std::ostream& out)
{
    typename traits<ValueType>::vector_type::HostMirror
            h_vec("h_vector",vec.extent(0));

    Kokkos::deep_copy(h_vec, vec);

    out << vec.label() << ": [";
    for (typename traits<ValueType>::vector_type::HostMirror::size_type i = 0;
         i < h_vec.extent(0) - 1;
         ++i) {
        out << h_vec(i) << ", ";
    }

    out << h_vec(h_vec.extent(0) - 1) << "]" << std::endl;
}









} // End namespace Vector

} // End namespace Impl

} // End namespace Numeric

} // End namespace ChipSum

#endif

#endif // __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
