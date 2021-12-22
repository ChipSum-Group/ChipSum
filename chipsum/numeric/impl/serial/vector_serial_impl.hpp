/*
 * @Description: 向量vector的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 16:11:28
 */

#ifndef __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_IMPL_HPP__



#include <fstream>

#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"




#include "vector_serial_dot_impl.hpp"
#include "vector_serial_axpby_impl.hpp"
#include "vector_serial_nrm1_impl.hpp"
#include "vector_serial_nrm2_impl.hpp"
#include "vector_serial_nrminf_impl.hpp"
#include "vector_serial_scal_impl.hpp"

/// 很抱歉的是vector的设计出现了失误，而且出现在比较棘手的kokkos后端上
/// 这部分工作我自己不会去做，主要的原因是暂时没有时间。
/// 这部分工作将交给后来的参与者。

/// 对于vector的具体修改内容介绍我会放在vector_kokkoskernels_impl.hpp里面

namespace ChipSum {
namespace Numeric {

template <typename ValueType, typename... Props>
struct Vector_Traits<ValueType,  ChipSum::Backend::Serial, Props...>
        : public Operator_Traits<ValueType> {
    using vector_type = typename ::std::vector<ValueType>;
    using size_type = typename ::std::vector<ValueType>::size_type;
    using scalar_type = ValueType;
    using value_type = ValueType;

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

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void create(::std::vector<ValueType> &dst,
                                    const ::std::size_t n) {
    dst.resize(n);
}

template <typename ST,typename ValueType>

CHIPSUM_FUNCTION_INLINE void create(::std::vector<ValueType> &dst,
                                    const ST& n,
                                    const ValueType *src
                                    ) {


    dst = ::std::vector<ValueType>(src, src + n);
}









template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void deep_copy(::std::vector<ValueType> &dst,
                                       const ::std::vector<ValueType> &src) {
    for(::std::size_t i=0;i<dst.size();++i){
        dst[i] = src[i];
    }
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void shallow_copy(::std::vector<ValueType> &dst,
                                          const ::std::vector<ValueType> &src) {
    dst = src;
}


template <typename ValueType,typename ST>

CHIPSUM_FUNCTION_INLINE ValueType &get_item(const ST& index,
                                            ::std::vector<ValueType> &vec) {


    return vec[index];
}

template <typename ValueType,typename OStreamT>

CHIPSUM_FUNCTION_INLINE void print(const ::std::vector<ValueType> &vec,
                                   OStreamT &out) {

    out << " [";
    for (::std::size_t i = 0; i < vec.size() - 1; ++i) {
        out << vec[i] << ", ";
    }

    out << vec[vec.size() - 1] << "]" << ::std::endl;
}


} // End namespace Vector
} // End namespace Impl

} // End namespace Numeric

} // End namespace ChipSum

#endif // End #ifndef HEADER_MACRO
