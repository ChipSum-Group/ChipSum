#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_SPMV_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_SPMV_IMPL_HPP__

#include "../../../chipsum_macro.h"

#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_default_types.hpp>

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Sparse {



template <typename ValueType,
          typename OrdinalType,
          typename SizeType>
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     const Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &b)
{
    KokkosSparse::spmv("N", 1, A, x,0, b);
}


template <typename ValueType,
          typename OrdinalType,
          typename SizeType>
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     const Kokkos::View<ValueType **> &B,
     Kokkos::View<ValueType **> &C) {
    KokkosSparse::spmv("N", static_cast<ValueType>(1), A, B,
                       static_cast<ValueType>(0), C);
}

template <typename ValueType,
          typename OrdinalType,
          typename SizeType,
          typename AlphaT,
          typename BetaT
          >
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &b,
     AlphaT alpha,
     BetaT beta) {
    KokkosSparse::spmv("N", alpha, A, x, beta, b);
}

template <typename ValueType,
          typename OrdinalType,
          typename SizeType,
          typename AlphaT,
          typename BetaT
          >
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     Kokkos::View<ValueType *> &x,
     Kokkos::View<ValueType *> &b,
     AlphaT alpha,
     BetaT beta,
     const char transA[]) {
    KokkosSparse::spmv(transA, alpha, A, x, beta, b);
}

}
}
}
}
#endif // CSR_KOKKOSKERNELS_SPMV_IMPL_HPP
