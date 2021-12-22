///
/// \file     vector_serial_scal_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-12-02
/// \brief    %stuff%
///
#ifndef __CHIPSUM_VECTOR_SERIAL_SCAL_IMPL_HPP__
#define __CHIPSUM_VECTOR_SERIAL_SCAL_IMPL_HPP__

#include <vector>
#include <algorithm>

#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Vector {



template <typename ValueType,typename AT>

CHIPSUM_FUNCTION_INLINE void scal(const ::std::vector<ValueType> &X,
                                  const AT& a,
                                  ::std::vector<ValueType> &R) {
    assert(R.size() == X.size());
    for (::std::size_t i = 0; i < X.size(); ++i) {
        R[i] = a * X[i];
    }

}



}
}
}
}

#endif // VECTOR_SERIAL_SCAL_IMPL_HPP
