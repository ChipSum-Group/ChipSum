#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include <type_traits>
#include "../backend/backend.hpp"

namespace ChipSum{
namespace Numeric{


template <typename ...Props>
struct Operator_Traits;

template <typename ScalarType>
struct Operator_Traits<ScalarType>{
    static_assert (std::is_floating_point<ScalarType>::value,"Float point error." );

    using nonconst_scalar_type = typename std::remove_const<ScalarType>::type;
    using const_scalar_type = typename
    std::enable_if<!std::is_const<ScalarType>::value,

    typename std::add_const<ScalarType>::type>::type;
};



template<typename ScalarType,typename SizeType>
struct Operator_Traits<ScalarType,SizeType>{

    static_assert (std::is_integral<SizeType>::value,"Integral type error." );

    using nonconst_scalar_type = typename Operator_Traits<ScalarType>::nonconst_scalar_type;

    using const_scalar_type = typename Operator_Traits<ScalarType>::const_scalar_type;

    using nonconst_size_type = typename std::remove_const<SizeType>::type;
    using const_size_type = typename
    std::enable_if<!std::is_const<SizeType>::value,

    typename std::add_const<SizeType>::type>::type;

};



template <typename ScalarType,typename SizeType,typename BackendType,typename ...Props >
struct Operator_Traits<ScalarType,SizeType,BackendType,Props...>{

    static_assert (std::is_integral<SizeType>::value,"Integral type error." );

    static_assert (std::is_base_of<ChipSum::Backend::Backend,
    BackendType>::value,"Parameter BackendType error." );



    using nonconst_scalar_type = typename Operator_Traits<ScalarType>::nonconst_scalar_type;

    using const_scalar_type = typename Operator_Traits<ScalarType>::const_scalar_type;

    using nonconst_size_type = typename
    Operator_Traits<ScalarType,SizeType>::nonconst_size_type;

    using const_size_type = typename
    Operator_Traits<ScalarType,SizeType>::const_size_type;

};





};
};

#endif // OPERATOR_HPP
