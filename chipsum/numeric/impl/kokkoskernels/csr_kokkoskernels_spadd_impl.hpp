#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_SPADD_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_SPADD_IMPL_HPP__

#include <string>

#include <KokkosSparse_spadd.hpp>
#include <KokkosKernels_default_types.hpp>

#include "../../../chipsum_macro.h"

namespace ChipSum {
namespace Numeric {
namespace Impl {
namespace Sparse {

template <typename ValueType,
          typename OrdinalType,
          typename SizeType>
CHIPSUM_FUNCTION_INLINE void
spadd(KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &C,
      const ValueType& alpha,
      KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
      const ValueType& beta,
      KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &B,
      const bool& isSortRows) {

    using device_type = typename Kokkos::Device<
    Kokkos::DefaultExecutionSpace,
    typename Kokkos::DefaultExecutionSpace::memory_space>;
    using execution_space = typename device_type::execution_space;
    using memory_space    = typename device_type::memory_space;

    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
    SizeType, OrdinalType, ValueType, execution_space, memory_space, memory_space>;

    KernelHandle kh;
    kh.create_spadd_handle(isSortRows);

    KokkosSparse::spadd_symbolic(&kh, A, B, C);
    KokkosSparse::spadd_numeric(&kh, alpha, A, beta, B, C);
    kh.destroy_spadd_handle();

    //KokkosKernels 列数输出还存在问题
    std::cout<<"Csr numCols: "<<C.numCols()<<std::endl;

}



}

}

}
}
#endif
