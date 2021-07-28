/* * * * * * * * * * * * * * * * * * * * *
*   File:     sparse_matrix.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__

#include <vector>
#include <type_traits>


#include "../chipsum_macro.h"
#include "numeric_traits.hpp"
#include "impl/crs_kokkoskernels_impl.hpp"
#include "impl/crs_serial_impl.hpp"
#include "vector.hpp"





namespace ChipSum {
namespace Numeric {

template<typename ...Props>
class SparseMatrix;

template<typename ScalarType,typename SizeType,typename SpFormat,typename BackendType,typename ...Props>
class SparseMatrix<ScalarType,SizeType,SpFormat,BackendType,Props...>{

public:

    using traits = Sparse_Traits<ScalarType,SizeType,SpFormat,BackendType,Props...>;
    using matrix_type = typename traits::matrix_format_type;
    using vector_type = Vector<ScalarType,SizeType,BackendType,Props...>;


private:

    matrix_type __data;


public:

    template<typename ...Args>
    CHIPSUM_FUNCTION_INLINE SparseMatrix(Args... args){
        ChipSum::Numeric::Impl::Sparse::Fill<ScalarType,SizeType,Props...>
                (__data,args...);
    }

    CHIPSUM_FUNCTION_INLINE vector_type
    operator*(const vector_type& v){
        typename vector_type::vector_type out_data;
        ChipSum::Numeric::Impl::Sparse::Spmv(__data,v,out_data);
        return rhs_type(out_data,v.GetSize());
    }


};

} // End namespace Numeric
} // End namespace ChipSum



#endif // SPARSEMATRIX_HPP
