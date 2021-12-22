#ifndef __CHIPSUM_VECTOR_SERIAL_NRM2_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_NRM2_IMPL_HPP__

#include <vector>
#include <cmath>

#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE ValueType norm2(const ::std::vector<ValueType> &X) {
    ValueType acc = 0.0;
    for (::std::size_t i = 0; i < X.size(); ++i) {
        acc += X[i] * X[i];
    }

    return ::std::sqrt(acc);
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void norm2(const ::std::vector<ValueType> &X,ValueType& acc) {
    acc = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        acc += X[i] *X[i];
    }

}

}
}
}
}

#endif // VECTOR_SERIAL_NRM2_IMPL_HPP
