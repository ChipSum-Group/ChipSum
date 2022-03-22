#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_SPTRSV_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_SPTRSV_IMPL_HPP__

#include <string>
#include <KokkosKernels_default_types.hpp>
#include "../../../chipsum_macro.h"
#include <Kokkos_DualView.hpp>
#include <KokkosSparse_sptrsv.hpp>

namespace ChipSum {
namespace Numeric {
namespace Impl {
namespace Sparse {

template <typename ValueType,
          typename OrdinalType,
          typename SizeType>
CHIPSUM_FUNCTION_INLINE void
sptrsv( KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
        Kokkos::DualView<ValueType *> &x,
        Kokkos::DualView<ValueType *> &b,
        bool &is_lower_tri) {
    using device_type = typename Kokkos::Device<
    Kokkos::DefaultExecutionSpace,
    typename Kokkos::DefaultExecutionSpace::memory_space>;
    using execution_space = typename device_type::execution_space;
    using memory_space    = typename device_type::memory_space;
    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
          SizeType, OrdinalType, ValueType, execution_space, memory_space, memory_space>;
    KernelHandle kh;
    using crsmat_t  = typename KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType>;

    const SizeType nrows = A.graph.numRows();

    // Create sptrsv kernel handle
    kh.create_sptrsv_handle(KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
    KokkosSparse::Experimental::sptrsv_symbolic(&kh, A.graph.row_map, A.graph.entries, A.values);
    Kokkos::fence();
    KokkosSparse::Experimental::sptrsv_solve(&kh, A.graph.row_map, A.graph.entries, A.values, b.d_view, x.d_view);
    Kokkos::fence();
    kh.destroy_sptrsv_handle();


}

}
}
}
}
#endif
