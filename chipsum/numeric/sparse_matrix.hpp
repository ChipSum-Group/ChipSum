/* * * * * * * * * * * * * * * * * * * * *
*   File:     sparse_matrix.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__


#include <type_traits>



#include "vector.hpp"

#include "impl/crs_serial_impl.hpp"
#include "impl/crs_kokkoskernels_impl.hpp"


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

    template<typename ...Args> // 用于应对不同格式稀疏矩阵的创建
    /**
     * @brief SparseMatrix
     * @param args
     */
    CHIPSUM_DECLARED_FUNCTION SparseMatrix(Args ...args){
        ChipSum::Numeric::Impl::Sparse::Create<ScalarType,SizeType>(__data,args...);
    }


    /**
     * @brief operator *
     * @param v
     * @return
     */
    CHIPSUM_FUNCTION_INLINE vector_type
    operator*(vector_type& v){
        vector_type ret(v.GetSize());
        ChipSum::Numeric::Impl::Sparse::Mult<ScalarType,SizeType>(
                    __data,v.GetData(),ret.GetData());
        return ret;
    }



};

} // End namespace Numeric
} // End namespace ChipSum



#endif // SPARSEMATRIX_HPP
