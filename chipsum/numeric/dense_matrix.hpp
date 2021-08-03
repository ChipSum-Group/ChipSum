#ifndef __CHIPSUM_DENSE_MATRIX_HPP__
#define __CHIPSUM_DENSE_MATRIX_HPP__

#include <fstream>

#include "numeric_traits.hpp"
#include "vector.hpp"
#include "impl/densemat_serial_impl.hpp"
#include "impl/densemat_kokkoskernels_impl.hpp"


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
    using vector_type = ChipSum::Numeric::Vector<ScalarType,SizeType,BackendType,Props...>;

private:

    matrix_type __data;
    size_type __nrow;
    size_type __ncol;

public:


    /**
     * @brief DenseMatrix
     * @param M
     * @param N
     */
    CHIPSUM_DECLARED_FUNCTION DenseMatrix(size_type M,size_type N)
        :__nrow(M),__ncol(N)
    {
        ChipSum::Numeric::Impl::DenseMat::Create<ScalarType,SizeType>(M,N,__data);
    }

    /**
     * @brief DenseMatrix
     * @param M
     * @param N
     * @param src
     */
    CHIPSUM_DECLARED_FUNCTION DenseMatrix(size_type M,size_type N,ScalarType* src)
        :__nrow(M),__ncol(N)
    {
        ChipSum::Numeric::Impl::DenseMat::Fill(M,N,src,__data);
    }

    /**
     * @brief GetData
     * @return
     */
    CHIPSUM_FUNCTION_INLINE matrix_type GetData(){return __data;}

    /**
     * @brief GetNRow
     * @return
     */
    CHIPSUM_FUNCTION_INLINE size_type GetRowNum(){return __nrow;}


    /**
     * @brief GetColNum
     * @return
     */
    CHIPSUM_FUNCTION_INLINE size_type GetColNum(){return __ncol;}


    /**
     * @brief operator *
     * @param v
     * @return
     */
    CHIPSUM_FUNCTION_INLINE matrix_type
    operator*(DenseMatrix& m){
        matrix_type ret;
        ChipSum::Numeric::Impl::DenseMat::Mult<ScalarType,SizeType>(__data,m.GetData(),ret);
        return ret;
    }

    template<typename ...Args>
    /**
     * @brief operator *
     * @param v
     * @return
     */
    CHIPSUM_FUNCTION_INLINE vector_type
    operator*(vector_type& v){
        vector_type ret(v.GetSize());
        ChipSum::Numeric::Impl::DenseMat::Mult<ScalarType,SizeType>(__nrow,__ncol,__data,v.GetData(),ret.GetData());
        return ret;
    }

    /**
     * @brief operator *
     * @param v
     * @return
     */
    CHIPSUM_FUNCTION_INLINE matrix_type
    operator*=(ScalarType s){
        ChipSum::Numeric::Impl::DenseMat::Scal<ScalarType,SizeType>(s,__data);
        return *this;
    }


    CHIPSUM_FUNCTION_INLINE ScalarType&
    operator()(const SizeType i,const SizeType j){
        return ChipSum::Numeric::Impl::DenseMat::GetItem<ScalarType,SizeType>(i,j,__nrow,__ncol,__data);
    }






    CHIPSUM_FUNCTION_INLINE void Print(std::ostream& out=std::cout){
        ChipSum::Numeric::Impl::DenseMat::Print<ScalarType,SizeType>(__nrow,__ncol,__data,out);
    }

};



} // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::DenseMatrix<double,std::size_t,ChipSum::Backend::KokkosKernels> Matrix;


#endif // __CHIPSUM_DENSE_MATRIX_HPP__
