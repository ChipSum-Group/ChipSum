///
/// \file     csr_serial_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-10
/// \brief    %stuff%
///

#ifndef __CHIPSUM_COO_SERIAL_IMPL_HPP__
#define __CHIPSUM_COO_SERIAL_IMPL_HPP__

#include <fstream>
#include <vector>
#include <cstdlib>
#include <assert.h>
#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"
#include "../../sparse_matrix_types.h"

// #include "../../../common/png_writer.hpp"
// #include "../../../common/bmp_writer.h"


namespace ChipSum {
namespace Numeric {


struct serial_coo_graph {
  ::std::vector<::std::size_t> row_map;
  ::std::vector<::std::size_t> col_map;
}__attribute__((aligned));

template <typename ValueType>
struct serial_coo_format {
  ::std::vector<ValueType> vals;
  serial_coo_graph graph;
  ::std::size_t row_num;
  ::std::size_t col_num;
  ::std::size_t nnz;
}__attribute__((aligned));


template <typename ValueType,typename ...Props>

struct Sparse_Traits<ValueType,
        ChipSum::Backend::Serial,
        ChipSum::Numeric::SparseTypes::Coo,
        Props...>
    : public Operator_Traits<ValueType> {

  using sp_type = serial_coo_format<ValueType>;
  using size_type = ::std::size_t;
  using ordinal_type = ::std::size_t;
  using backend_type = ChipSum::Backend::Serial;
  using format_type = ChipSum::Numeric::SparseTypes::Coo;

  using graph_type = serial_coo_graph;
  using row_map_type = ::std::vector<::std::size_t>;
  using col_map_type = ::std::vector<::std::size_t>;
  using values_type = ::std::vector<ValueType>;
  using value_type = ValueType;
};

namespace Impl {

namespace Sparse {

template <typename ValueType, typename S1, typename S2>
// 通过POD数据创建CSR格式矩阵
CHIPSUM_FUNCTION_INLINE void
create(serial_coo_format<ValueType> &A,
       S1 nrows,
       S1 ncols,
       S1 annz,
       S2 *row_map,
       S2 *col_map,
       ValueType *values) {

  //标准COO数据排序
  for(S1 i = 0; i < annz; ++i){
    for(S1 j = 0; j < annz - i - 1; ++j){
      if(row_map[j] > row_map[j+1]){
        std::swap(row_map[j],row_map[j+1]);
        std::swap(col_map[j],col_map[j+1]);
        std::swap(values[j],values[j+1]);
      }else if(row_map[j] == row_map[j+1]){
        if(col_map[j] > col_map[j+1]){
          std::swap(col_map[j],col_map[j+1]);
          std::swap(values[j],values[j+1]);
        }
      }
    }
  }

  A.vals = std::vector<ValueType>(values, values + annz);
  A.graph.row_map = std::vector<::std::size_t>(row_map, row_map + annz);
  A.graph.col_map = std::vector<::std::size_t>(col_map, col_map + annz);
  A.row_num = nrows;
  A.col_num = ncols;
  A.nnz = annz;
}

template <typename ValueType, typename S>
// 创建未初始化的COO格式矩阵
CHIPSUM_FUNCTION_INLINE void create(serial_coo_format<ValueType> &A,
                                    S nrows,
                                    S ncols,
                                    S annz
                                    ) {

  A.vals.resize(annz);
  A.graph.row_map.resize(annz);
  A.graph.col_map.resize(annz);
  A.row_num = nrows;
  A.col_num = ncols;
  A.nnz = annz;
}

// COO矩阵插入元素
template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void insert(serial_coo_format<ValueType> &A,
                                    std::size_t row,
                                    std::size_t col,
                                    ValueType val){
  assert(row<A.row_num);
  assert(col<A.col_num);
  size_t off_row = std::lower_bound(A.graph.row_map.begin(), A.graph.row_map.end(), row) - A.graph.row_map.begin();
  size_t row_elems = std::upper_bound(A.graph.row_map.begin(), A.graph.row_map.end(), row) -
                     A.graph.row_map.begin() - off_row;
  if(row_elems != 0){
    //检查插入的数据是否存在
    auto iter_begin = A.graph.col_map.begin() + off_row;
    auto iter_end = iter_begin + row_elems;
    auto iter = std::find(iter_begin, iter_end, col);
    if(iter != iter_end){
      std::cerr<<"the inserted data already exists"<<std::endl;
      abort();
    }

    size_t off_val = std::lower_bound(iter_begin, iter_end, col) -
                     A.graph.col_map.begin();
    A.graph.row_map.insert(A.graph.row_map.begin() + off_val, row);
    A.graph.col_map.insert(A.graph.col_map.begin() + off_val, col);
    A.vals.insert(A.vals.begin() + off_val, val);
  } else{
    A.graph.row_map.insert(A.graph.row_map.begin() + off_row, row);
    A.graph.col_map.insert(A.graph.col_map.begin() + off_row, col);
    A.vals.insert(A.vals.begin() + off_row, val);
  }
  ++A.nnz;
}

template <typename ValueType, typename S>
// 由COO创建CSR
CHIPSUM_FUNCTION_INLINE void get_csr_data(
                                    serial_coo_format<ValueType> &A,
                                    std::vector<S> &csr_row_map,
                                    std::vector<S> &csr_col_map,
                                    std::vector<ValueType> &values
                                    ) {
  csr_row_map.resize(A.row_num+1);
  csr_col_map.resize(A.nnz);
  values.resize(A.nnz);

  for (size_t i = 0; i < A.nnz; ++i)
    {
        values[i] = A.vals[i];
        csr_col_map[i] = A.graph.col_map[i];
        csr_row_map[A.graph.row_map[i] + 1]++;
    }

  for (size_t i = 1; i < A.row_num + 1; i++)
  {
      csr_row_map[i] += csr_row_map[i-1];
  }

}


} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
