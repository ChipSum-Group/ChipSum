#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__


#include <KokkosSparse.hpp>
#include <KokkosSparse_spmv.hpp>

#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"
#include "../../chipsum_macro.h"
#include

namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
struct Sparse_Traits<ScalarType,SizeType,Csr,ChipSum::Backend::KokkosKernels,Props...>
        : public Operator_Traits<ScalarType,SizeType>{


    using matrix_format_type = KokkosSparse::CrsMatrix<ScalarType,SizeType>;


};



CHIPSUM_FUNCTION_INLINE Spmv()


}
}



#endif // CRS_KOKKOSKERNELS_IMPL_H
