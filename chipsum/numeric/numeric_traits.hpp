#ifndef __NUMERIC_TRAITS_HPP__
#define __NUMERIC_TRAITS_HPP__

#include <type_traits>
#include "../backend/backend.hpp"

namespace ChipSum{
namespace Numeric{


template <typename ...Props>
struct Operator_Traits;



//template<typename ScalarType,typename SizeType>
//struct Operator_Traits<ScalarType,SizeType>{



//};



template <typename ScalarType,typename SizeType,typename BackendType,typename ...Props >
struct Operator_Traits<ScalarType,SizeType,BackendType,Props...>{

    static_assert (std::is_scalar<ScalarType>::value,"[ERR] scalar type error" );

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
    BackendType>::value,"Parameter BackendType error." );




};


template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
struct Vector_Traits: public Operator_Traits<ScalarType,SizeType,BackendType>
{
    using vector_type = void;
    using vector_type_reference = void;
    using const_vector_type_reference = void;

//    using traits1 = Operator_Traits<ScalarType,Props...>;

//    using nonconst_scalar_type = typename traits1::nonconst_scalar_type;
//    using const_scalar_type = typename traits1::const_scalar_type;
//    using nonconst_scalar_type_reference = typename traits1::nonconst_scalar_type_reference;
//    using const_scalar_type_reference = typename traits1::const_scalar_type_reference;


//    using traits2 = Operator_Traits<ScalarType,SizeType,Props...>;
//    using nonconst_size_type = typename traits2::nonconst_size_type;
//    using const_size_type = typename traits2::const_size_type;
//    using nonconst_size_type_reference = typename traits2::nonconst_size_type_reference;
};


template<typename ScalarType,typename SizeType,typename ...Props>
struct Sparse_Traits: public Operator_Traits<ScalarType,SizeType,Props...>
{

    using matrix_format_type = void;

    using traits1 = Operator_Traits<ScalarType,Props...>;

    using nonconst_scalar_type = typename traits1::nonconst_scalar_type;
    using const_scalar_type = typename traits1::const_scalar_type;
    using nonconst_scalar_type_reference = typename traits1::nonconst_scalar_type_reference;
    using const_scalar_type_reference = typename traits1::const_scalar_type_reference;


    using traits2 = Operator_Traits<ScalarType,SizeType,Props...>;
    using nonconst_size_type = typename traits2::nonconst_size_type;
    using const_size_type = typename traits2::const_size_type;
    using nonconst_size_type_reference = typename traits2::nonconst_size_type_reference;
};



}
}

#endif // OPERATOR_HPP
