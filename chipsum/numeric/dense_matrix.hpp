///
/// \file     dense_matrix.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-01
/// \brief    稠密矩阵用户接口
///

#ifndef __CHIPSUM_DENSE_MATRIX_HPP__
#define __CHIPSUM_DENSE_MATRIX_HPP__

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/densemat_kokkoskernels_impl.hpp"
#endif
#include "impl/densemat_serial_impl.hpp"
#include "numeric_traits.hpp"
#include "scalar.hpp"
#include "vector.hpp"





namespace ChipSum {
namespace Numeric {


template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
class DenseMatrix<ScalarType, SizeType, BackendType, Props...> {

public:
  using traits =
      DenseMatrix_Traits<ScalarType, SizeType, BackendType, Props...>;
  using matrix_type = typename traits::matrix_type;
  using matrix_type_reference =
      typename std::add_lvalue_reference<matrix_type>::type;
  using const_matrix_type_reference =
      typename std::add_const<matrix_type_reference>::type;

  using size_type = typename traits::size_type;
  using const_size_type = typename traits::const_size_type;

  using vector_type =
      ChipSum::Numeric::Vector<ScalarType, SizeType, BackendType, Props...>;

private:
  matrix_type __data;
  size_type __nrow;
  size_type __ncol;

public:

  ///
  /// \brief DenseMatrix 构造一个M行N列的稠密矩阵，该矩阵未初始化
  /// \param M M 行数
  /// \param N N 列数
  ///
  CHIPSUM_DECLARED_FUNCTION DenseMatrix(const_size_type M, const_size_type N)
      : __nrow(M), __ncol(N) {
    ChipSum::Numeric::Impl::DenseMat::create<ScalarType, SizeType>(M, N,
                                                                   __data);
  }


  ///
  /// \brief DenseMatrix 构造一个M行N列的稠密矩阵，并用src赋值该矩阵
  /// \param M M 行数
  /// \param N N 列数
  /// \param src 原数据
  ///
  CHIPSUM_DECLARED_FUNCTION DenseMatrix(const_size_type M, const_size_type N,
                                        ScalarType *src)
      : __nrow(M), __ncol(N) {
    ChipSum::Numeric::Impl::DenseMat::create<ScalarType, SizeType>(M, N,
                                                                   __data);
    ChipSum::Numeric::Impl::DenseMat::fill<ScalarType, SizeType>(M, N, src,
                                                                 __data);

  }


  ///
  /// \brief GetData 获取矩阵数据
  /// \return 后端数据
  ///
  CHIPSUM_FUNCTION_INLINE const_matrix_type_reference GetData() {
    return __data;
  }


  ///
  /// \brief GetRowNum 获取矩阵行数
  /// \return  矩阵行数
  ///
  CHIPSUM_FUNCTION_INLINE size_type GetRowNum() { return __nrow; }


  ///
  /// \brief GetColNum 获取矩阵列数
  /// \return 矩阵列数
  ///
  CHIPSUM_FUNCTION_INLINE size_type GetColNum() { return __ncol; }


  ///
  /// \brief operator * GEMM
  /// \param m 稠密矩阵
  /// \return 稠密矩阵（结果）
  ///
  CHIPSUM_FUNCTION_INLINE DenseMatrix operator*(DenseMatrix &m) {
    DenseMatrix ret(__nrow, m.GetColNum());
    ChipSum::Numeric::Impl::DenseMat::mult<ScalarType, SizeType>(
        __nrow,m.GetColNum(),m.GetRowNum(), __data,m.GetData(), ret.GetData());
    return ret;
  }

  ///
  /// \brief GEMM C=A*B 当C为已初始化的矩阵时，强烈建议采用此接口进行GEMM运算
  /// \param B  参与运算的另一矩阵
  /// \param C  结果
  ///
  CHIPSUM_FUNCTION_INLINE void GEMM(DenseMatrix &B,DenseMatrix& C) {
    ChipSum::Numeric::Impl::DenseMat::mult<ScalarType, SizeType>(
        __nrow,B.GetColNum(),B.GetRowNum(), __data,B.GetData(), C.GetData());

  }

  template <typename... Args>
  ///
  /// \brief operator * GEMV
  /// \param v 向量
  /// \return 向量（结果）
  ///
  CHIPSUM_FUNCTION_INLINE vector_type operator*(vector_type &v) {
    vector_type ret(__ncol);
    ChipSum::Numeric::Impl::DenseMat::mult<ScalarType, SizeType>(
        __nrow, __ncol, __data, v.GetData(), ret.GetData());
    return ret;
  }


  ///
  /// \brief operator *= A*=a
  /// \param a 系数
  /// \return A（结果）
  ///
  CHIPSUM_FUNCTION_INLINE matrix_type operator*=(const ScalarType a) {
    ChipSum::Numeric::Impl::DenseMat::scal<ScalarType, SizeType>(a, __data);
    return *this;
  }


  ///
  /// \brief operator () 获取A(i,j)
  /// \param i 行索引
  /// \param j 列索引
  /// \return A(i,j)
  ///
  CHIPSUM_FUNCTION_INLINE ScalarType &operator()(const_size_type i,
                                                 const_size_type j) {
    return ChipSum::Numeric::Impl::DenseMat::get_item<ScalarType, SizeType>(
        i, j, __nrow, __ncol, __data);
  }



  ///
  /// \brief operator () 获取A(i,j)（只读）
  /// \param i 行数
  /// \param j 列数
  /// \return A(i,j)
  ///
  CHIPSUM_FUNCTION_INLINE const ScalarType &operator()(const_size_type i,
                                                 const_size_type j) const{
    return ChipSum::Numeric::Impl::DenseMat::get_item<ScalarType, SizeType>(
        i, j, __nrow, __ncol, __data);
  }


  ///
  /// \brief Print 打印函数
  /// \param out 输出流
  ///
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::DenseMat::print<ScalarType, SizeType>(
        __nrow, __ncol, __data, out);
  }
};

} // End namespace Numeric
} // End namespace ChipSum

///
/// \brief Matrix 默认的
///
typedef ChipSum::Numeric::DenseMatrix<CSFloat,
                                      CSInt,
                                      ChipSum::Backend::DefaultBackend>
    Matrix;

#endif // __CHIPSUM_DENSE_MATRIX_HPP__
