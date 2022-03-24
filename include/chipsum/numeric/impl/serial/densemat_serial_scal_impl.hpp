#ifndef __CHIPSUM_DENSEMAT_SERIAL_SCAL_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_SERIAL_SCAL_IMPL_HPP__

#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {


template <typename ValueType,typename DT>
// A = alpha*A
CHIPSUM_FUNCTION_INLINE void scal(const DT &A,
                                  DT &B,
                                  const ValueType& alpha)
{
    for (std::size_t i = 0; i < A.data.size(); ++i) {
        B.data[i] = alpha*A.data[i];
    }
}


}
}
}
}



#endif // DENSEMAT_SERIAL_SCAL_IMPL_HPP
