#ifndef __CHIPSUM_VECTOR_SERIAL_AXPBY_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_AXPBY_IMPL_HPP__

#include <vector>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void add(
        const ::std::vector<ValueType> &X,
        ::std::vector<ValueType> &Y) {
    assert(X.size() == Y.size());
    for (::std::size_t i = 0; i < Y.size(); ++i) {
        Y[i] +=  X[i];
    }
}

template <typename ValueType,typename AlphaT>

CHIPSUM_FUNCTION_INLINE void axpby(
        const ::std::vector<ValueType> &X,
        ::std::vector<ValueType> &Y,
        const AlphaT& a) {
    assert(X.size() == Y.size());
    for (::std::size_t i = 0; i < Y.size(); ++i) {
        Y[i] += a * X[i];
    }
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void axpby(const ::std::vector<ValueType> &X,
                                   ::std::vector<ValueType> &Y,
                                   const AlphaT& a,
                                   const BetaT& b) {
    assert(X.size() == Y.size());

    for (::std::size_t i = 0; i < Y.size(); ++i) {
        Y[i] = a * X[i] + b * Y[i];
    }
}
}}}}
#endif // VECTOR_SERIAL_AXPBY_IMPL_HPP
