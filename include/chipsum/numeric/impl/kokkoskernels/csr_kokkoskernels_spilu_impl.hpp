#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_SPILU_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_SPILU_IMPL_HPP__

#include <string>
#include <KokkosSparse_spiluk.hpp>
#include <KokkosKernels_default_types.hpp>
#include "../../../chipsum_macro.h"
#define EXPAND_FACT 6 // a factor used in expected sizes of L and U

namespace ChipSum {
namespace Numeric {
namespace Impl {
namespace Sparse {

template <typename ValueType,
          typename OrdinalType,
          typename SizeType>
CHIPSUM_FUNCTION_INLINE void
spilu( KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &A,
       KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &L,
       KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType> &U,
       SizeType &fill_lev) {
    using device_type = typename Kokkos::Device<
    Kokkos::DefaultExecutionSpace,
    typename Kokkos::DefaultExecutionSpace::memory_space>;
    using execution_space = typename device_type::execution_space;
    using memory_space    = typename device_type::memory_space;
    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
          SizeType, OrdinalType, ValueType, execution_space, memory_space, memory_space>;
    KernelHandle kh;
    using crsmat_t  = typename KokkosSparse::CrsMatrix<ValueType,OrdinalType, default_device,void,SizeType>;
    using graph_t         = typename crsmat_t::StaticCrsGraphType;
    using lno_view_t      = typename graph_t::row_map_type::non_const_type;
    using lno_nnz_view_t  = typename graph_t::entries_type::non_const_type;
    using scalar_view_t   = typename crsmat_t::values_type::non_const_type;

    const SizeType N = A.graph.numRows();
    const SizeType nnzA = A.graph.entries.extent(0);
    kh.create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_RP, N,
                            EXPAND_FACT*nnzA*(fill_lev+1), EXPAND_FACT*nnzA*(fill_lev+1));
    auto spiluk_handle = kh.get_spiluk_handle();
    
    lno_view_t     L_row_map("L_row_map", N + 1);
    lno_nnz_view_t L_entries("L_entries", spiluk_handle->get_nnzL());
    scalar_view_t  L_values ("L_values",  spiluk_handle->get_nnzL());
    lno_view_t     U_row_map("U_row_map", N + 1);
    lno_nnz_view_t U_entries("U_entries", spiluk_handle->get_nnzU());
    scalar_view_t  U_values ("U_values",  spiluk_handle->get_nnzU());	

    KokkosSparse::Experimental::spiluk_symbolic(&kh, fill_lev, A.graph.row_map, A.graph.entries, 
                                                L_row_map, L_entries, U_row_map, U_entries);

    Kokkos::resize(L_entries, spiluk_handle->get_nnzL());
    Kokkos::resize(L_values,  spiluk_handle->get_nnzL());
    Kokkos::resize(U_entries, spiluk_handle->get_nnzU());
    Kokkos::resize(U_values,  spiluk_handle->get_nnzU());	
    
    KokkosSparse::Experimental::spiluk_numeric(&kh, fill_lev, A.graph.row_map, A.graph.entries, A.values, 
                                               L_row_map, L_entries, L_values, U_row_map, U_entries, U_values);
    L=crsmat_t("L", N, N, spiluk_handle->get_nnzL(), L_values, L_row_map, L_entries);
    U=crsmat_t("U", N, N, spiluk_handle->get_nnzU(), U_values, U_row_map, U_entries);

    kh.destroy_spiluk_handle();


}

}
}
}
}
#endif
