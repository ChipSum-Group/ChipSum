/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 16:03:13
 */


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

template <typename... Props> class DenseMatrix;

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
  /**
   * @description: 构造一个M行N列的稠密矩阵，该矩阵未初始化
   * @param {size_type} M 行数
   * @param {size_type} N 列数
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION DenseMatrix(const_size_type M, const_size_type N)
      : __nrow(M), __ncol(N) {
    ChipSum::Numeric::Impl::DenseMat::create<ScalarType, SizeType>(M, N,
                                                                   __data);
  }

  /**
   * @description: 构造一个M行N列的稠密矩阵，并用src赋值该矩阵
   * @param {size_type} M 行数
   * @param {size_type} N 列数
   * @param {ScalarType*} src 原数据
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION DenseMatrix(const_size_type M, const_size_type N,
                                        ScalarType *src)
      : __nrow(M), __ncol(N) {
    ChipSum::Numeric::Impl::DenseMat::create<ScalarType, SizeType>(M, N,
                                                                   __data);
    ChipSum::Numeric::Impl::DenseMat::fill<ScalarType, SizeType>(M, N, src,
                                                                 __data);
                                                                
  }

  /**
   * @description: 获取矩阵数据
   * @param {*}
   * @return {const_matrix_type_reference} 稠密矩阵数据
   */
  CHIPSUM_FUNCTION_INLINE const_matrix_type_reference GetData() {
    return __data;
  }

  /**
   * @description: 获取矩阵行数
   * @param {*}
   * @return {size_type} 矩阵行数
   */
  CHIPSUM_FUNCTION_INLINE size_type GetRowNum() { return __nrow; }

  /**
   * @description: 获取矩阵列数
   * @param {*} 
   * @return {*} 矩阵列数
   */
  CHIPSUM_FUNCTION_INLINE size_type GetColNum() { return __ncol; }

  /**
   * @description: GEMM
   * @param {DenseMatrix} m 稠密矩阵 
   * @return {*} 稠密矩阵（结果）
   */
  CHIPSUM_FUNCTION_INLINE DenseMatrix operator*(DenseMatrix &m) {
    DenseMatrix ret(__nrow, m.GetColNum());
    ChipSum::Numeric::Impl::DenseMat::mult<ScalarType, SizeType>(
        __nrow,m.GetColNum(),m.GetRowNum(), __data,m.GetData(), ret.GetData());
    return ret;
  }

  template <typename... Args>
  /**
   * @description: GEMV
   * @param {vector_type} v 向量
   * @return {*} 向量（结果）
   */
  CHIPSUM_FUNCTION_INLINE vector_type operator*(vector_type &v) {
    vector_type ret(__ncol);
    ChipSum::Numeric::Impl::DenseMat::mult<ScalarType, SizeType>(
        __nrow, __ncol, __data, v.GetData(), ret.GetData());
    return ret;
  }

  /**
   * @description: A*=a
   * @param {const ScalarType} a 系数
   * @return {*} A（结果）
   */
  CHIPSUM_FUNCTION_INLINE matrix_type operator*=(const ScalarType a) {
    ChipSum::Numeric::Impl::DenseMat::scal<ScalarType, SizeType>(a, __data);
    return *this;
  }
  /**
   * @description: 获取Aij
   * @param {const_size_type} i 行索引
   * @param {const_size_type} j 列索引
   * @return {*} Aij
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE ScalarType &operator()(const_size_type i,
                                                 const_size_type j) {
    return ChipSum::Numeric::Impl::DenseMat::get_item<ScalarType, SizeType>(
        i, j, __nrow, __ncol, __data);
  }

    /**
   * @description: 获取Aij（只读）
   * @param {*} i 行数
   * @param {*} j 列数
   * @return {*} Aij
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE const ScalarType &operator()(const_size_type i,
                                                 const_size_type j) const{
    return ChipSum::Numeric::Impl::DenseMat::get_item<ScalarType, SizeType>(
        i, j, __nrow, __ncol, __data);
  }

  /**
   * @description: 打印函数
   * @param {std::ostream&} 输出流
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::DenseMat::print<ScalarType, SizeType>(
        __nrow, __ncol, __data, out);
  }
};

} // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::DenseMatrix<double, std::size_t,
                                      ChipSum::Backend::DefaultBackend>
    Matrix;

#endif // __CHIPSUM_DENSE_MATRIX_HPP__
