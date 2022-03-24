#ifndef __CHIPSUM_VECTOR_SERIAL_DOT_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_DOT_IMPL_HPP__

#include <vector>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void dot(const ::std::vector<ValueType> &x,

                                 const ::std::vector<ValueType> &y,

                                 ValueType &r) {

    assert(x.size() == y.size());

    for (::std::size_t i = 0; i < x.size(); ++i) {
        r += x[i] * y[i];
    }
}


template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void dot(const ::std::vector<ValueType> &x,
                                 const ::std::vector<ValueType> &y) {

    ValueType r = 0;
    assert(x.size() == y.size());

    for (::std::size_t i = 0; i < x.size(); ++i) {
        r += x[i] * y[i];
    }

    return r;
}


}
}
}
}

#endif // VECTOR_SERIAL_DOT_IMPL_HPP
