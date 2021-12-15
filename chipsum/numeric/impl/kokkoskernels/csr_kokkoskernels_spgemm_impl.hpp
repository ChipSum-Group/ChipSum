#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_SPGEMM_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_SPGEMM_IMPL_HPP__

#include <string>

#include <KokkosSparse_spgemm.hpp>
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
spgemm(KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
       KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &B,
       KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &C) {
    using device_type = typename Kokkos::Device<
    Kokkos::DefaultExecutionSpace,
    typename Kokkos::DefaultExecutionSpace::memory_space>;
    using execution_space = typename device_type::execution_space;
    using memory_space    = typename device_type::memory_space;

    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
    SizeType, OrdinalType, ValueType, execution_space, memory_space, memory_space>;
    KernelHandle kh;
    kh.set_team_work_size(16);
    kh.set_dynamic_scheduling(true);

    // Select an spgemm algorithm, limited by configuration at compile-time and
    // set via the handle Some options: {SPGEMM_KK_MEMORY, SPGEMM_KK_SPEED,
    // SPGEMM_KK_MEMSPEED, /*SPGEMM_CUSPARSE, */ SPGEMM_MKL}
    ::std::string myalg("SPGEMM_KK_MEMORY");
    KokkosSparse::SPGEMMAlgorithm spgemm_algorithm =
            KokkosSparse::StringToSPGEMMAlgorithm(myalg);
    kh.create_spgemm_handle(spgemm_algorithm);

    KokkosSparse::spgemm_symbolic(kh, A, false, B, false, C);
    KokkosSparse::spgemm_numeric(kh, A, false, B, false, C);



}

}

}

}
}
#endif
