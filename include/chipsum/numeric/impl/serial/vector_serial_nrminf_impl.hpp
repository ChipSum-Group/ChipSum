#ifndef __CHIPSUM_VECTOR_SERIAL_NRMINF_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_NRMINF_IMPL_HPP__
#include <vector>
#include <algorithm>

#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {





template <typename ValueType>

CHIPSUM_FUNCTION_INLINE ValueType norminf(const ::std::vector<ValueType> &X) {


    return (*::std::max_element(X.begin(),X.end()));
}


template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void norminf(const ::std::vector<ValueType> &X,ValueType& r) {


    r = (*::std::max_element(X.begin(),X.end()));

}

}
}
}
}
#endif // VECTOR_SERIAL_NRMINF_IMPL_HPP
