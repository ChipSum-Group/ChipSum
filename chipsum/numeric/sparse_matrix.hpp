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
#include "dense_matrix.hpp"

#include "impl/crs_kokkoskernels_impl.hpp"
#include "impl/crs_serial_impl.hpp"




namespace ChipSum {
namespace Numeric {


template<typename ...Props>
class SparseMatrix;

enum GSAlgorithm{DEFAULT, PERMUTED, TEAM, CLUSTER, TWOSTAGE};

template<typename ScalarType,typename SizeType,typename SpFormat,typename BackendType,typename ...Props>
class SparseMatrix<ScalarType,SizeType,SpFormat,BackendType,Props...>{

public:

    using traits = Sparse_Traits<ScalarType,SizeType,SpFormat,BackendType,Props...>;
    using sp_type = typename traits::sp_type;
    using size_type = typename traits::size_type;


    using vector_type = Vector<ScalarType,SizeType,BackendType,Props...>;
    using dense_type = DenseMatrix<ScalarType,SizeType,BackendType,Props...>;


private:

    sp_type __data;
    size_type __nrow;
    size_type __ncol;
    size_type __annz;




public:

    template<typename ...Args> // 用于应对不同格式稀疏矩阵的创建
    /**
     * @brief SparseMatrix
     * @param args
     */
    CHIPSUM_DECLARED_FUNCTION SparseMatrix(size_type nrow,
                                           size_type ncol,
                                           size_type annz,
                                           Args ...args)
        :__nrow(nrow),__ncol(ncol),__annz(annz)
    {
        ChipSum::Numeric::Impl::Sparse::Create<ScalarType,SizeType>(nrow,ncol,annz,__data,args...);
    }


    /**
     * @brief GetColNum
     * @return
     */
    CHIPSUM_FUNCTION_INLINE size_type GetColNum(){return __ncol;}


    /**
     * @brief GetColNum
     * @return
     */
    CHIPSUM_FUNCTION_INLINE size_type GetRowNum(){return __nrow;}

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


    /**
     * @brief operator *
     * @param m
     * @return
     */
    CHIPSUM_FUNCTION_INLINE dense_type
    operator*(dense_type& m){
        dense_type ret(__nrow,m.GetColNum());
        ChipSum::Numeric::Impl::Sparse::Mult<ScalarType,SizeType>(
                    __data,m.GetData(),ret.GetData());
        return ret;
    }



    CHIPSUM_FUNCTION_INLINE void
    GSSmooth(GSAlgorithm algo = DEFAULT)
    {
        ChipSum::Numeric::Impl::Sparse::GaussSeidelSmooth
                <ScalarType,SizeType>
                (__data,__nrow,__ncol,static_cast<char>(algo));
    }

};

} // End namespace Numeric
} // End namespace ChipSum



#endif // SPARSEMATRIX_HPP
