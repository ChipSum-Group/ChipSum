#ifndef __CHIPSUM_DENSE_MATRIX_HPP__
#define __CHIPSUM_DENSE_MATRIX_HPP__

#include <type_traits>



#include "numeric_traits.hpp"
#include "vector.hpp"

#include "impl/densemat_blas_impl.hpp"


namespace ChipSum {
namespace Numeric {

template<typename ...Props>
class DenseMatrix;

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
class DenseMatrix<ScalarType,SizeType,BackendType,Props...>{

public:

    using traits = DenseMatrix_Traits<ScalarType,SizeType,BackendType,Props...>;
    using matrix_type = typename traits::matrix_type;
    using size_type = typename traits::size_type;
    using vector_type = Vector<ScalarType,SizeType,BackendType,Props...>;


private:

    matrix_type __data;


public:

    template<typename ...Args> // 用于应对不同格式稀疏矩阵的创建
    /**
     * @brief SparseMatrix
     * @param args
     */
    CHIPSUM_FUNCTION_INLINE DenseMatrix(Args ...args){
        // TODO
    }


    CHIPSUM_FUNCTION_INLINE matrix_type GetData(){
        return __data;
    }

    /**
     * @brief operator *
     * @param v
     * @return
     */
    CHIPSUM_FUNCTION_INLINE vector_type
    operator*(DenseMatrix& m){
//        ChipSum::Numeric::Impl::DenseMat::Mult(__data,m.GetData());
    }



};

} // End namespace Numeric
} // End namespace ChipSum




#endif // __CHIPSUM_DENSE_MATRIX_HPP__
