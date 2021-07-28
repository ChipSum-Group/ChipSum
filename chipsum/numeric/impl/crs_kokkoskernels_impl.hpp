/* * * * * * * * * * * * * * * * * * * * *
*   File:     crs_kokkoskernels_impl.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__


#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_trsv.hpp>

#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"
#include "../../chipsum_macro.h"


static int spm_name = 0;

namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Sparse_Traits<ScalarType,SizeType,Csr,ChipSum::Backend::KokkosKernels,Props...>
        : public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::KokkosKernels,Props...>{


    using matrix_format_type = KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>;

    using graph_type = typename matrix_format_type::staticcrsgraph_type;
    using row_map_type = typename matrix_format_type::row_map_type;
    using col_map_type = typename matrix_format_type::index_type;
    using matrix_values_type = typename matrix_format_type::values_type;


};


namespace Impl {

namespace Sparse {



template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Fill：//TODO
 * @param nrows
 * @param ncols
 * @param annz
 * @param row_map
 * @param col_map
 * @param values
 * @param A
 */
CHIPSUM_FUNCTION_INLINE void Fill(
        KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& A,
        const SizeType nrows,
        const SizeType ncols,
        const SizeType annz,
        SizeType* row_map,
        SizeType* col_map,
        ScalarType* values

        )
{

    using crs_t = KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>;

    A.ctor_impl("spm_"+std::to_string(spm_name),
                static_cast<typename crs_t::ordinal_type>(nrows),
                static_cast<typename crs_t::ordinal_type>(ncols),
                static_cast<typename crs_t::size_type>(annz),
                static_cast<typename crs_t::value_type*>(values),
                static_cast<typename crs_t::ordinal_type*>(row_map),
                static_cast<typename crs_t::ordinal_type*>(col_map));

}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Spmv
 * @param A
 * @param x
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void Spmv(
        KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& A,
        Kokkos::View<ScalarType*>& x,
        Kokkos::View<ScalarType*>& b)
{
    KokkosSparse::spmv("N",static_cast<ScalarType>(1.0),A,x,static_cast<ScalarType>(0.0),b);
}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Spmv
 * @param alpha
 * @param A
 * @param x
 * @param beta
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void Spmv(
        ScalarType alpha,
        KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& A,
        Kokkos::View<ScalarType*>& x,
        ScalarType beta,
        Kokkos::View<ScalarType*> b)
{
    KokkosSparse::spmv("N",alpha,A,x,beta,b);
}


} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
