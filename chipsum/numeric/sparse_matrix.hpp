/*
 * @Description: 稀疏矩阵用户接口
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 15:18:14
 */

#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__

#include <type_traits>
#include <fstream>

#include "dense_matrix.hpp"
#include "vector.hpp"

#include "impl/crs_serial_impl.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/crs_kokkoskernels_impl.hpp"
#endif


namespace ChipSum {
namespace Numeric {

template <typename... Props> class SparseMatrix;

enum GSAlgorithm { DEFAULT, PERMUTED, TEAM, CLUSTER, TWOSTAGE };

template <typename ScalarType, typename SizeType, typename SpFormat,
          typename BackendType, typename... Props>
class SparseMatrix<ScalarType, SizeType, SpFormat, BackendType, Props...> {

public:
  using traits =
      Sparse_Traits<ScalarType, SizeType, SpFormat, BackendType, Props...>;
  using sp_type = typename traits::sp_type;
  using size_type = typename traits::size_type;

  using vector_type = Vector<ScalarType, SizeType, BackendType, Props...>;
  using dense_type = DenseMatrix<ScalarType, SizeType, BackendType, Props...>;

private:
  sp_type __data;
  size_type __nrow;
  size_type __ncol;
  size_type __annz;

public:
  template <typename... Args> // 用于应对不同格式稀疏矩阵的创建
  /**
   * @description:
   * @param {size_type} nrow 行数
   * @param {size_type} ncol 列数
   * @param {size_type} annz 非零元数
   * @param {Args} 稀疏矩阵所需的其他数据
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION SparseMatrix(size_type nrow, size_type ncol,
                                         size_type annz, Args... args)
      : __nrow(nrow), __ncol(ncol), __annz(annz) {
    ChipSum::Numeric::Impl::Sparse::create<ScalarType, SizeType>(
        nrow, ncol, annz, __data, args...);
  }
  

  /**
   * @description:
   * @param {*}
   * @return {size_type} 获取列数
   */
  CHIPSUM_FUNCTION_INLINE size_type GetColNum() { return __ncol; }

  /**
   * @description:
   * @param {*}
   * @return {size_type} 获取行数
   */
  CHIPSUM_FUNCTION_INLINE size_type GetRowNum() { return __nrow; }

  /**
   * @description: SpMV
   * @param {vector_type} 向量
   * @return {*} 向量（结果）
   */
  CHIPSUM_FUNCTION_INLINE vector_type operator*(vector_type &v) {
    vector_type ret(v.GetSize());
    ChipSum::Numeric::Impl::Sparse::mult<ScalarType, SizeType>(
        __nrow,__ncol,__data, v.GetData(), ret.GetData());
    return ret;
  }

  /**
   * @description: SpMD（稀疏矩阵乘稠密矩阵）
   * @param {dense_type} 稠密矩阵
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE dense_type operator*(dense_type &m) {
    dense_type ret(__nrow, m.GetColNum());
    ChipSum::Numeric::Impl::Sparse::mult<ScalarType, SizeType>(
        __nrow,m.GetColNum(),__ncol,__data, m.GetData(), ret.GetData());
    return ret;
  }


  /**
   * @description: 打印（调试用）
   * @param {std::ostream&} 输出流
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream& out=std::cout)
  {
    ChipSum::Numeric::Impl::Sparse::print<ScalarType,SizeType>(__data,out);
  }



  /**
   * @description: 打印pattern（调试用）
   * @param {std::ostream&} 输出流
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE void PrintPattern(std::ostream& out=std::cout)
  {
    ChipSum::Numeric::Impl::Sparse::print_pattern<ScalarType,SizeType>(__data,out);
  }


  
  /**
   * @description: 保存pattern为图片（调试用）
   * @param {std::ostream&} 输出流
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE void SavePatternFig(const char* filename)
  {
    // 还有一些BUG:BMP格式的稀疏矩阵必须是方阵，且行列数必须是4的倍数。（待修复）
    ChipSum::Numeric::Impl::Sparse::save_figure<ScalarType,SizeType>(__data,filename);
  }

};

} // End namespace Numeric
} // End namespace ChipSum
typedef ChipSum::Numeric::SparseMatrix<double, size_t,
                                       ChipSum::Numeric::SparseTypes::Csr,
                                       ChipSum::Backend::DefaultBackend>
    CSR;
#endif // SPARSEMATRIX_HPP
