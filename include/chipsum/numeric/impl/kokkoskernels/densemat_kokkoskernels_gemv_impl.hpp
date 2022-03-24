#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMV_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_GEMV_IMPL_HPP__

#include <KokkosBlas2_gemv.hpp>
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
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y)
{

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*n, hipblas A is n*m
     // x is n, hipblas x is n
     // Thus A*x is m*1. Hipblas A_T*x is m*1, and then transformer to m*1
     int m = A.extent(0);
     int n = A.extent(1);

     int lda = n;
     int incx = 1;
     int incy = 1;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemv(handle, rocblas_operation_transpose, n, m, 1.0, A.d_view.data(), lda, 
                    x.d_view.data(), incx, 0.0, y.d_view.data(), incy);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemv("N", 1, A.d_view, x.d_view,
                     0, y.d_view);
#endif
}


template <typename ValueType>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y,
     const char transA[])
{

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*n, hipblas A is n*m
     // x is n, hipblas x is n
     // Thus A*x is m*1. Hipblas A_T*x is m*1, and then transformer to m*1
     int m = A.extent(0);
     int n = A.extent(1);

     int lda = n;
     int incx = 1;
     int incy = 1;

     rocblas_operation transa = (transA=='N') ? rocblas_operation_transpose : rocblas_operation_none;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemv(handle, transa, n, m, 1.0, A.d_view.data(), lda, 
                    x.d_view.data(), incx, 0.0, y.d_view.data(), incy);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemv(transA, 1, A.d_view, x.d_view,
                     0, y.d_view);
#endif
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y,
     const AlphaT& a,
     const BetaT& b
     )
{

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*n, hipblas A is n*m
     // x is n, hipblas x is n
     // Thus A*x is m*1. Hipblas A_T*x is m*1, and then transformer to m*1
     int m = A.extent(0);
     int n = A.extent(1);

     int lda = n;
     int incx = 1;
     int incy = 1;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemv(handle, rocblas_operation_transpose, n, m, a, A.d_view.data(), lda, 
                    x.d_view.data(), incx, b, y.d_view.data(), incy);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemv("N", a, A.d_view, x.d_view,
                     b, y.d_view);
#endif
}

template <typename ValueType,typename AlphaT,typename BetaT>

CHIPSUM_FUNCTION_INLINE void
gemv(const Kokkos::DualView<ValueType **> &A,
     const Kokkos::DualView<ValueType *> &x,
     Kokkos::DualView<ValueType *> &y,
     const AlphaT& a,
     const BetaT& b,
     const char transA[])
{

#ifdef CHIPSUM_USE_HIPBLAS
     // hipblas is col major
     // A is m*n, hipblas A is n*m
     // x is n, hipblas x is n
     // Thus A*x is m*1. Hipblas A_T*x is m*1, and then transformer to m*1
     int m = A.extent(0);
     int n = A.extent(1);

     int lda = n;
     int incx = 1;
     int incy = 1;

     rocblas_operation transa = (transA=='N') ? rocblas_operation_transpose : rocblas_operation_none;

     rocblas_handle handle;
     rocblas_create_handle(&handle);
     rocblas_dgemv(handle, transa, n, m, a, A.d_view.data(), lda, 
                    x.d_view.data(), incx, b, y.d_view.data(), incy);
     rocblas_destroy_handle(handle);
#else
     KokkosBlas::gemv(transA, a, A.d_view, x.d_view,
                     b, y.d_view);
#endif
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_GEMM_IMPL_HPP
