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
using namespace std;
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

/// 由于直接用view<type*>作为vector的数据底层，导致了数据无法灵活在
/// device和host之间流转。所以这部分的底层数据类型替换为dualview。


namespace ChipSum {
namespace Numeric {

template <typename ValueType,typename ...Props>
struct Vector_Traits<ValueType, ChipSum::Backend::KokkosKernels,Props...>
        : public Operator_Traits<ValueType> {
    using vector_type = Kokkos::DualView<ValueType *>;
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
        Kokkos::DualView<ValueType *>& x,
        const ST& n
        )
{

    x.resize(n);
}

template <typename ValueType,typename ST>
CHIPSUM_FUNCTION_INLINE void create(
        Kokkos::DualView<ValueType *>& x,
        const ST& n,
        ValueType* src
        ) {
    
    Kokkos::resize(x.d_view,n);
    x.h_view = typename Kokkos::View<ValueType*>::HostMirror(src,n);
    Kokkos::deep_copy(x.d_view, x.h_view);
 
}

template <typename ValueType,typename ST>
CHIPSUM_FUNCTION_INLINE void create(
        Kokkos::DualView<ValueType *>& x,
        const ST& n,
        ValueType value
        ) {
    x.resize(n);
    Kokkos::deep_copy(x.h_view, value);
    Kokkos::deep_copy(x.d_view, x.h_view);
}



template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
deep_copy(
        const Kokkos::DualView<ValueType *>& dst,
        const Kokkos::DualView<ValueType *>& src)
{

    Kokkos::deep_copy(dst.h_view, src.h_view);
    Kokkos::deep_copy(dst.d_view, src.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
device_to_host(Kokkos::DualView<ValueType *>& x) 
{
    Kokkos::deep_copy(x.h_view, x.d_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
host_to_device(Kokkos::DualView<ValueType *>& x)
{ 
    Kokkos::deep_copy(x.d_view, x.h_view);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE
typename ::std::add_lvalue_reference<
typename traits<ValueType>::vector_type::value_type
>::type
get_item(
        const Kokkos::DualView<ValueType *>& vec,
        const size_t index
        )
{
    return vec.h_view(index);
}

template <typename ValueType>
CHIPSUM_SPECIAL_INLINE
ValueType&
item(
        const Kokkos::DualView<ValueType *>& vec,
        const size_t index
        )
{
    return vec.d_view(index);
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE
void print( Kokkos::DualView<ValueType *>& vec,
           ::std::ostream& out)
{
    Kokkos::deep_copy(vec.h_view, vec.d_view);
    out<< vec.h_view.label()<<": ";
    for(typename traits<ValueType>::vector_type::size_type i = 0; i < vec.extent(0)-1; ++i)
        {
            out << vec.h_view(i) << ", ";
        }
        out<<vec.h_view(vec.extent(0)-1)<<endl;

}









} // End namespace Vector

} // End namespace Impl

} // End namespace Numeric

} // End namespace ChipSum

#endif

#endif // __CHIPSUM_VECTOR_KOKKOSKERNELS_IMPL_HPP__
