/*
 * @Description: 标量scalar的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 13:58:38
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-16 10:28:17
 */


#ifndef __CHIPSUM_SCALAR_SERIAL_IMPL_HPP__
#define __CHIPSUM_SCALAR_SERIAL_IMPL_HPP__

#include <fstream>
#include <vector>

#include "../../../backend/backend.hpp"
#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"

namespace ChipSum {
namespace Numeric {

template <typename ValueType,typename... Props>
struct Scalar_Traits<ValueType,ChipSum::Backend::Serial, Props...>
        : public Operator_Traits<ValueType> {
    using scalar_type = ValueType;
    using value_type = ValueType;
    using host_type = ValueType;
};

namespace Impl {
namespace Scalar {
template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void create(ValueType &) {
    return;
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void create(const ValueType s, ValueType &r) {
    r = s;
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void deep_copy(ValueType &r,const ValueType& s) {
    r = s;
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE ValueType& get_item(ValueType& s) {
    return s;
}

template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void get_item(const ValueType s, ValueType &r) {
    r = s;
}


template <typename ValueType,typename OStreamT>

CHIPSUM_FUNCTION_INLINE void print(const ValueType s, OStreamT &out) {

    cout << "scalar_serial: " << s << endl;
}
} // namespace Scalar
} // namespace Impl

} // namespace Numeric
} // namespace ChipSum

#endif
