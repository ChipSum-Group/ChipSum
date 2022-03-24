#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_SPMV_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_SPMV_IMPL_HPP__

#include "../../../chipsum_macro.h"

#include <KokkosSparse_spmv.hpp>
#include <Kokkos_DualView.hpp>
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
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &b)
{
    KokkosSparse::spmv("N", 1, A, x.d_view,0, b.d_view);
}


template <typename ValueType,
          typename OrdinalType,
          typename SizeType>
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C) {
    KokkosSparse::spmv("N", static_cast<ValueType>(1), A, B.d_view,
                       static_cast<ValueType>(0), C.d_view);
}

template <typename ValueType,
          typename OrdinalType,
          typename SizeType,
          typename AlphaT,
          typename BetaT
          >
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &b,
     AlphaT alpha,
     BetaT beta) {
    KokkosSparse::spmv("N", alpha, A, x.d_view, beta, b.d_view);
}

template <typename ValueType,
          typename OrdinalType,
          typename SizeType,
          typename AlphaT,
          typename BetaT
          >
CHIPSUM_FUNCTION_INLINE void
spmv(const KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
     Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &b,
     AlphaT alpha,
     BetaT beta,
     const char transA[]) {
    KokkosSparse::spmv(transA, alpha, A, x.d_view, beta, b.d_view);
}

}
}
}
}
#endif // CSR_KOKKOSKERNELS_SPMV_IMPL_HPP
