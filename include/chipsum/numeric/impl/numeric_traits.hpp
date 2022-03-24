
///
/// \file     numeric_traits.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-23
/// \brief    %stuff%
///

#ifndef __CHIPSUM_NUMERIC_TRAITS_HPP__
#define __CHIPSUM_NUMERIC_TRAITS_HPP__

#include <type_traits>
#include "../backend/backend.hpp"

namespace ChipSum{
namespace Numeric{


template <typename ...Props>
struct Operator_Traits;


template <typename ScalarType,typename ...Props >
struct Operator_Traits<ScalarType,Props...>{
    static_assert (std::is_scalar<ScalarType>::value,"[ERR] template parameter ScalarType error" );
};




template<typename ...Props>
struct Scalar_Traits: public Operator_Traits<Props...>
{};


template<typename ...Props>
struct Vector_Traits: public Operator_Traits<Props...>
{};

template<typename ...Props>
struct DenseMatrix_Traits: public Operator_Traits<Props...>
{};


template<typename ...Props>
struct Sparse_Traits: public Operator_Traits<Props...>
{};

template<size_t NDIM, typename ...Props>
struct Tensor_Traits: public Operator_Traits<Props...>
{};


//template <typename... Props>
//class Scalar;

//template<typename ...Props>
//class Vector;

//template<typename ...Props>
//class DenseMatrix;

//template<typename ...Props>
//class SparseMatrix;

}
}

#endif // OPERATOR_HPP
