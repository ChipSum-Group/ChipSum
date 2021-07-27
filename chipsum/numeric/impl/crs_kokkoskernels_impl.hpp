#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__


#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse.hpp>
#include <KokkosSparse_spmv.hpp>

#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"
#include "../../chipsum_macro.h"


namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Sparse_Traits<ScalarType,SizeType,Csr,ChipSum::Backend::KokkosKernels,Props...>
        : public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::KokkosKernels,Props...>{


    using matrix_format_type = KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>;


};


template<typename ScalarType,typename SizeType>
CHIPSUM_FUNCTION_INLINE void Spmv(KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& m,){}


}
}



#endif // CRS_KOKKOSKERNELS_IMPL_H
