#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__


#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse.hpp>
#include <KokkosSparse_spmv.hpp>

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


};


namespace Impl {

namespace Sparse {



template<typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void Fill(
        const SizeType nrows,
        const SizeType ncols,
        const SizeType annz,
        SizeType* row_map,
        SizeType* col_map,
        ScalarType* values,
        KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& A
        )
{


    A.ctor_impl("spm_"+std::to_string(spm_name),nrows,ncols,annz,values,row_map,col_map);

}


template<typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void Spmv(
        KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& A,
        Kokkos::View<ScalarType*>& x,
        Kokkos::View<ScalarType*>& b)
{
    KokkosSparse::spmv("N",static_cast<ScalarType>(1.0),A,x,static_cast<ScalarType>(0.0),b);
}


template<typename ScalarType,typename SizeType,typename ...Props>
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
