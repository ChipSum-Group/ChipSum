#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP__

#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_DualView.hpp>

#include "../../../chipsum_macro.h"

#ifdef CHIPSUM_USE_HIPBLAS
#include <hip/hip_runtime.h>
#include <rocblas.h>
#endif


namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace DenseMat {

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C) {

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*k, hipblas A is k*m
     // B is k*n, hipblas B is n*k
     // Thus A*B is m*n. Hipblas B*A is n*m, and then transformer to m*n
     int m = A.extent(0);
     int n = B.extent(1);
     int k = A.extent(1);

     int lda = n;
     int ldb = k;
     int ldc = n;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_none, n, m, k, 1.0, B.d_view.data(), lda, 
                    A.d_view.data(), ldb, 0.0, C.d_view.data(), ldc);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemm("N", "N", 1, A.d_view, B.d_view,
                         0, C.d_view);
#endif
}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C,
     const char transA[],
     const char transB[]
     ) {

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*k, hipblas A is k*m
     // B is k*n, hipblas B is n*k
     // Thus A*B is m*n. Hipblas B*A is n*m, and then transformer to m*n
     int m = A.extent(0);
     int n = B.extent(1);
     int k = A.extent(1);

     int lda = n;
     int ldb = k;
     int ldc = n;

     rocblas_operation transa = (transA=='N') ? rocblas_operation_none : rocblas_operation_transpose;
     rocblas_operation transb = (transB=='N') ? rocblas_operation_none : rocblas_operation_transpose;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemm(handle, transb, transa, n, m, k, 1.0, B.d_view.data(), lda, 
                    A.d_view.data(), ldb, 0.0, C.d_view.data(), ldc);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemm(transA, transB, 1, A.d_view, B.d_view,
                     0, C.d_view);
#endif
}

template <typename ValueType,typename AlphaT,typename BetaT>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C,
     const AlphaT& a,
     const BetaT& b
     ) {

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*k, hipblas A is k*m
     // B is k*n, hipblas B is n*k
     // Thus A*B is m*n. Hipblas B*A is n*m, and then transformer to m*n
     int m = A.extent(0);
     int n = B.extent(1);
     int k = A.extent(1);

     int lda = n;
     int ldb = k;
     int ldc = n;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_none, n, m, k, a, B.d_view.data(), lda, 
                    A.d_view.data(), ldb, b, C.d_view.data(), ldc);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemm("N", "N", a, A.d_view, B.d_view,
                     b, C.d_view);
#endif
}

template <typename ValueType,typename AlphaT,typename BetaT>
CHIPSUM_FUNCTION_INLINE void
gemm(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType **> &B,
     Kokkos::DualView<ValueType **> &C,
     const AlphaT& a,
     const BetaT& b,
     const char transA[],
     const char transB[]
     ) {

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*k, hipblas A is k*m
     // B is k*n, hipblas B is n*k
     // Thus A*B is m*n. Hipblas B*A is n*m, and then transformer to m*n
     int m = A.extent(0);
     int n = B.extent(1);
     int k = A.extent(1);

     int lda = n;
     int ldb = k;
     int ldc = n;

     rocblas_operation transa = (transA=='N') ? rocblas_operation_none : rocblas_operation_transpose;
     rocblas_operation transb = (transB=='N') ? rocblas_operation_none : rocblas_operation_transpose;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemm(handle, transb, transa, n, m, k, a, B.d_view.data(), lda, 
                    A.d_view.data(), ldb, b, C.d_view.data(), ldc);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemm(transA, transB, a, A.d_view, B.d_view,
                     b, C.d_view);
#endif
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP
