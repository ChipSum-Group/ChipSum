#ifndef __CHIPSUM_VECTOR_SERIAL_NRM1_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_NRM1_IMPL_HPP__

#include <vector>

#include "../../../chipsum_macro.h"



namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE ValueType norm1(const ::std::vector<ValueType> &X) {
    ValueType acc = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        acc += X[i];
    }
    return acc;
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void norm1(const ::std::vector<ValueType> &X,ValueType& acc) {
    acc = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        acc += X[i];
    }

}

}
}
}
}

#endif // VECTOR_SERIAL_NRM1_IMPL_HPP
