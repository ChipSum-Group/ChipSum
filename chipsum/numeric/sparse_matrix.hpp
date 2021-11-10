///
/// \file     sparse_matrix.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-01
/// \brief    稀疏矩阵用户接口
///

#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__

#include <type_traits>
#include <fstream>

#include "dense_matrix.hpp"
#include "vector.hpp"

#include "impl/csr_serial_impl.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/csr_kokkoskernels_impl.hpp"
#endif


namespace ChipSum {
namespace Numeric {




///
/// \brief 稀疏矩阵用户接口
///
/// \tparam ScalarType 元素的数据类型，一般为浮点数或者整数。
///         暂时不支持内置的Scalar和complex类型
///
/// \tparam OrdinalType 行、列所用的整数描述类型。一般来说单
///         节点的行数和列数不会比行邻接表大，常见情况是小很多。
///         所以在这里我参考Trilinos对他们的整型描述符进行了
///         参数分离
///
/// \tparam SizeType Row Map在可能会很大，甚至比行、列大很多
///         。（参考Trilinos）
///
/// \tparam SpFormat 稀疏矩阵的类型。可以参考sparse_matrix_types.h
///
///
/// \tparam BackendType 后端类型。参考ChipSum User Manual
///
template <typename ScalarType, typename OrdinalType,
          typename SizeType, typename SpFormat,
          typename BackendType, typename... Props>

class SparseMatrix<ScalarType, OrdinalType,SizeType,
        SpFormat, BackendType, Props...>
{

public:
    using traits =
    Sparse_Traits<ScalarType, OrdinalType, SizeType, SpFormat, BackendType, Props...>;
    using sp_type = typename traits::sp_type;
    using size_type = typename traits::size_type;

    using ordinal_type = OrdinalType;

    using vector_type = Vector<ScalarType, SizeType, BackendType, Props...>;
    using dense_type = DenseMatrix<ScalarType, SizeType, BackendType, Props...>;

private:

    sp_type __data;
    ordinal_type __nrow;
    ordinal_type __ncol;
    size_type __annz;

public:

    SparseMatrix()=default;

    template <typename... Args> // 用于应对不同格式稀疏矩阵的创建

    ///
    /// \brief SparseMatrix 稀疏矩阵模板。SpFormat表示稀疏矩阵的储存格式，至少会支持CSC/CSR/COO，后续会添加一些新的后端。
    /// \param nrow 行数
    /// \param ncol 列数
    /// \param annz 非零元数
    /// \param args 稀疏矩阵所需的其他数据
    ///
    CHIPSUM_DECLARED_FUNCTION SparseMatrix(ordinal_type nrow, ordinal_type ncol,
                                           size_type annz, Args... args)
        : __nrow(nrow), __ncol(ncol), __annz(annz) {
        ChipSum::Numeric::Impl::Sparse::create<ScalarType,OrdinalType, SizeType>(
                    nrow, ncol, annz, __data, args...);
    }


    ///
    /// \brief GetData 获取后端类型数据
    /// \return 后端类型数据
    ///
    CHIPSUM_FUNCTION_INLINE sp_type& GetData() {return __data;}


    ///
    /// \brief GetColNum 获取列数
    /// \return 列数
    ///
    CHIPSUM_FUNCTION_INLINE ordinal_type GetColNum() { return __ncol; }

    ///
    /// \brief GetRowNum 获取行数
    /// \return 行数
    ///
    CHIPSUM_FUNCTION_INLINE ordinal_type GetRowNum() { return __nrow; }


    ///
    /// \brief GetNNZ 获取非零元数
    /// \return 非零元数
    ///
    CHIPSUM_FUNCTION_INLINE size_type& GetNNZ() {return __annz;}

    ///
    /// \brief operator * SpMV operator*版SpMV主要是
    ///                   方便用户创造结果向量，性能并不
    ///                   好
    /// \param x 向量
    /// \return 向量（结果）
    ///
    CHIPSUM_FUNCTION_INLINE vector_type operator*(vector_type &x) {
        vector_type ret(x.GetSize());
        ChipSum::Numeric::Impl::Sparse::mult<ScalarType,OrdinalType, SizeType>(
                    __nrow,__ncol,__data, x.GetData(), ret.GetData());
        return ret;
    }


    ///
    /// \brief operator * SpMV（稀疏矩阵乘稠密矩阵）
    /// \param m 稠密矩阵
    /// \return
    ///
    CHIPSUM_FUNCTION_INLINE dense_type operator*(dense_type &m) {
        dense_type ret(__nrow, m.GetColNum());
        ChipSum::Numeric::Impl::Sparse::mult<ScalarType, OrdinalType,SizeType>(
                    __nrow,m.GetColNum(),__ncol,__data, m.GetData(), ret.GetData());
        return ret;
    }




    //  ///
    //  /// \brief operator * SPGEMM（稀疏矩阵乘稀疏矩阵）
    //  /// \param m
    //  /// \return
    //  ///
    //  CHIPSUM_FUNCTION_INLINE SparseMatrix operator*(SparseMatrix& m) {

    //      return nullptr;
    //  }


    ///
    /// \brief SpMV 在需要考虑性能的时候，强烈建议采用此接口进行稀疏矩阵乘向量
    /// \param x 左端项
    /// \param y 右端项
    ///
    CHIPSUM_FUNCTION_INLINE void Multiply(vector_type &x,vector_type &y)
    {
        ChipSum::Numeric::Impl::Sparse::mult<ScalarType,OrdinalType, SizeType>(
                    __nrow,__ncol,__data, x.GetData(), y.GetData());
    }

    ///
    /// \brief SpMV 在需要考虑性能的时候，强烈建议采用此接口进行稀疏矩阵乘矩阵
    /// \param x 左端项
    /// \param y 右端项
    ///
//    CHIPSUM_FUNCTION_INLINE void Multiply(dense_type &x,dense_type &y)
//    {
//        ChipSum::Numeric::Impl::Sparse::mult<ScalarType, OrdinalType,SizeType>(
//                    __nrow,__ncol,__data, x.GetData(), y.GetData());
//    }

    ///
    /// \brief Print 打印（调试用）
    /// \param out 输出流
    ///
    CHIPSUM_FUNCTION_INLINE void Print(std::ostream& out=std::cout)
    {
        ChipSum::Numeric::Impl::Sparse::print<ScalarType,SizeType>(__data,out);
    }


    ///
    /// \brief PrintPattern 打印pattern（调试用）
    /// \param out 输出流
    ///
    CHIPSUM_FUNCTION_INLINE void PrintPattern(std::ostream& out=std::cout)
    {
        ChipSum::Numeric::Impl::Sparse::print_pattern<ScalarType,SizeType>(__data,out);
    }


    ///
    /// \brief SavePatternFig 保存pattern为图片（调试用）
    /// \param filename 文件名
    ///
    CHIPSUM_FUNCTION_INLINE void SavePatternFig(const char* filename)
    {
        // 还有一些BUG:BMP格式的稀疏矩阵必须是方阵，且行列数必须是4的倍数。（待修复）
        ChipSum::Numeric::Impl::Sparse::save_figure<ScalarType,SizeType>(__data,filename);
    }


    ///
    /// \brief Multiply SPGEMM C=A*B 算子有待测试
    /// \param B [IN]
    /// \param C [OUT]
    ///
    CHIPSUM_FUNCTION_INLINE void Multiply(SparseMatrix& B,SparseMatrix& C){
        ChipSum::Numeric::Impl::Sparse::spgemm<ScalarType,OrdinalType,SizeType>
                (__data,B.GetData(),C.GetData());

    }

};

} // End namespace Numeric
} // End namespace ChipSum
typedef
ChipSum::Numeric::SparseMatrix<CSFloat, CSInt, CSInt,
ChipSum::Numeric::SparseTypes::Csr,
ChipSum::Backend::DefaultBackend>
CSR;
#endif // SPARSEMATRIX_HPP
