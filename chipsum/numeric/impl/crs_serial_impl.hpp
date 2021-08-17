/*
 * @Description: CRS矩阵的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-17 15:54:57
 */

#ifndef __CHIPSUM_CRS_SERIAL_IMPL_HPP__
#define __CHIPSUM_CRS_SERIAL_IMPL_HPP__

#include <fstream>
#include <vector>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"

namespace ChipSum {
namespace Numeric {

template <typename SizeType, typename... Props> struct StaticGraph {
  ::std::vector<SizeType> row_map;
  ::std::vector<SizeType> col_map;
};

template <typename ScalarType, typename SizeType, typename... Props>
struct CrsFormat {
  ::std::vector<ScalarType> vals;
  StaticGraph<SizeType> graph;
  ::std::size_t col_num;
};

template <typename ScalarType, typename SizeType, typename... Props>

struct Sparse_Traits<ScalarType, SizeType, SparseTypes::Csr,
                     ChipSum::Backend::Serial, Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Serial,
                             Props...> {

  using sp_type = CrsFormat<ScalarType, SizeType>;
  using size_type = ::std::size_t;

  using graph_type = StaticGraph<SizeType>;
  using row_map_type = ::std::vector<SizeType>;
  using col_map_type = ::std::vector<SizeType>;
  using values_type = ::std::vector<ScalarType>;
};

namespace Impl {

namespace Sparse {

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 创建初始化的CRS矩阵
 * @param {*} nrows 行数
 * @param {*} ncols 列数
 * @param {*} annz 非零元数
 * @param {*} A 稀疏矩阵（out）
 * @param {*} row_map 行邻接表
 * @param {*} col_map 列邻接表
 * @param {*} values 非零元
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Create(const SizeType nrows, const SizeType ncols, const SizeType annz,
       CrsFormat<ScalarType, SizeType> &A, SizeType *row_map, SizeType *col_map,
       ScalarType *values) {
  CHIPSUM_UNUSED(ncols);
  A.vals = ::std::vector<ScalarType>(values, values + annz);
  A.graph.row_map = ::std::vector<SizeType>(row_map, row_map + nrows + 1);
  A.graph.col_map = ::std::vector<SizeType>(col_map, col_map + annz);
  A.col_num = ncols;
}

template <typename ScalarType, typename SizeType, typename... Props>
/** ------UNFINISHED------
 * @description: 创建未初始化的CRS矩阵
 * @param {*} A 稀疏矩阵（out）
 * @param {*} row_map_size 行邻接表长度
 * @param {*} col_map_size 列邻接表长度
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Create(CrsFormat<ScalarType, SizeType> &A,
                                    const ::std::size_t row_map_size,
                                    const ::std::size_t col_map_size) {
  CHIPSUM_UNUSED(row_map_size);
  A.vals.resize(col_map_size);
  A.graph.row_map.resize(col_map_size);
  A.graph.col_map.resize(col_map_size);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: b=Ax（SpMV）
 * @param {*} A 稀疏矩阵
 * @param {*} x 向量
 * @param {*} b 向量（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Mult(::std::size_t M, ::std::size_t N, CrsFormat<ScalarType, SizeType> &A,
     ::std::vector<ScalarType> &x, ::std::vector<ScalarType> &b) {

  assert(M == b.size());

  for (::std::size_t i = 0; i < M; ++i)
    b[i] = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = A.graph.row_map[i];
    ::std::size_t end = A.graph.row_map[i + 1];
    for (::std::size_t j = 0; j < end - start; ++j) {
      b[i] += A.vals[start + j] * x[A.graph.col_map[start + j]];
    }
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: b = alpha*Ax+beta*b
 * @param {*} alpha 稀疏矩阵A的系数
 * @param {*} A 稀疏矩阵
 * @param {*} x 向量
 * @param {*} beta 向量b的系数
 * @param {*} b 向量（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Mult(ScalarType alpha, CrsFormat<ScalarType, SizeType> &A,
     ::std::vector<ScalarType> &x, ScalarType beta, ::std::vector<ScalarType> &b) {

  for (::std::size_t i = 0; i < b.size(); ++i) {
    ::std::size_t start = A.graph.row_map[i];
    ::std::size_t end = A.graph.row_map[i + 1];
    for (::std::size_t j = 0; j < end - start; ++j) {
      b[i] += beta * b[i] +
              alpha * A.vals[start + j] * x[A.graph.col_map[start + j]];
    }
  }
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印出稀疏矩阵的数据信息，多用于调试
 * @param {*} A 稀疏矩阵
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void Print(CrsFormat<ScalarType, SizeType> &A,
                                   ::std::ostream &out) {
  out << "spm_serial:"
      << "(rows=" << A.graph.row_map.size() - 1 << ", entries=" << A.vals.size()
      << ")" << ::std::endl;

  out << "rowmap: [";
  for (::std::size_t i = 0; i < A.graph.row_map.size() - 1; ++i) {
    out << A.graph.row_map[i] << ",";
  }
  out << A.graph.row_map[A.graph.row_map.size() - 1] << "]" << ::std::endl;

  out << "entries: [";
  for (::std::size_t i = 0; i < A.graph.col_map.size() - 1; ++i) {
    out << A.graph.col_map[i] << ",";
  }
  out << A.graph.col_map[A.graph.col_map.size() - 1] << "]" << ::std::endl;

  out << "values: [";
  for (::std::size_t i = 0; i < A.vals.size() - 1; ++i) {
    out << A.vals[i] << ",";
  }
  out << A.vals[A.vals.size() - 1] << "]" << ::std::endl;
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: 打印出稀疏矩阵的Pattern信息，多用于调试
 * @param {*} A 稀疏矩阵
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void PrintPattern(CrsFormat<ScalarType, SizeType> &A,
                                   ::std::ostream &out) {
  
  ::std::size_t M = A.graph.row_map.size()-1;
  ::std::size_t N = A.col_num;

  ::std::size_t entry_cnt = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = A.graph.row_map[i];
    for (::std::size_t j = 0; j < N; ++j) {
      char info = 'o';
      if(A.graph.col_map[start+entry_cnt]==j){
         info = '+';
         ++entry_cnt;
      }
      out<<info<<" ";
    }
    entry_cnt = 0;
    out<<::std::endl;
  }
}


} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
