/*
 * @Description: 数值类数据的共性提取和特性分离
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:34:55
 */



#ifndef __CHIPSUM_NUMERIC_TRAITS_HPP__
#define __CHIPSUM_NUMERIC_TRAITS_HPP__

#include <type_traits>
#include "../backend/backend.hpp"

namespace ChipSum{
namespace Numeric{


template <typename ...Props>
struct Operator_Traits;


template <typename ScalarType,typename SizeType,typename BackendType,typename ...Props >
struct Operator_Traits<ScalarType,SizeType,BackendType,Props...>{

    static_assert (std::is_scalar<ScalarType>::value,"[ERR] template parameter ScalarType error" );
    static_assert (std::is_integral<SizeType>::value,"[ERR] template parameter SizeType error" );

    using nonconst_scalar_type = typename std::remove_const<ScalarType>::type;

    using const_scalar_type = typename
    std::enable_if<!std::is_const<ScalarType>::value,
    typename std::add_const<ScalarType>::type>::type;

    using nonconst_scalar_type_reference = typename std::add_lvalue_reference_t<nonconst_scalar_type>;
    using const_scalar_type_reference = typename std::add_lvalue_reference_t<const_scalar_type>;

    using nonconst_scalar_type_pointer = typename std::add_pointer_t<nonconst_scalar_type>;
    using const_scalar_type_pointer = typename std::add_const_t<nonconst_scalar_type_pointer>;


    static_assert (std::is_integral<SizeType>::value,"[ERR] integral type error" );

    using nonconst_size_type = typename std::remove_const<SizeType>::type;
    using const_size_type = typename
    std::enable_if<!std::is_const<SizeType>::value,
    typename std::add_const<SizeType>::type>::type;


    using nonconst_size_type_reference = typename std::add_lvalue_reference_t<nonconst_size_type>;
    using const_size_type_reference = typename std::add_const_t<nonconst_size_type_reference>;

    using nonconst_size_type_pointer = typename std::add_pointer_t<nonconst_size_type>;
    using const_size_type_pointer = typename std::add_const_t<nonconst_size_type_pointer>;

    static_assert (std::is_base_of<ChipSum::Backend::BackendBase,
    BackendType>::value,"[ERR] Parameter BackendType error." );




};

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
struct Scalar_Traits: public Operator_Traits<ScalarType,SizeType,BackendType,Props...>
{
    using scalar_type = void;

    using size_type = void;

};


template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
struct Vector_Traits: public Operator_Traits<ScalarType,SizeType,BackendType,Props...>
{
    using vector_type = void;

    using size_type = void;

};

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
struct DenseMatrix_Traits: public Operator_Traits<ScalarType,SizeType,BackendType,Props...>
{
    using matrix_type = void;

    using size_type = void;

};


template<typename ScalarType,typename SizeType,typename SparseType,typename BackendType,typename ...Props>
struct Sparse_Traits: public Operator_Traits<ScalarType,SizeType,BackendType,Props...>
{

    using matrix_format_type = void;
    using graph_type = void;

    using matrix_values_type = void;


};




template <typename... Props>
class Scalar;

template<typename ...Props>
class Vector;

template<typename ...Props>
class DenseMatrix;

template<typename ...Props>
class SparseMatrix;

}
}

#endif // OPERATOR_HPP
