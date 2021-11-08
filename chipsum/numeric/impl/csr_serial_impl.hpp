/*
 * @Description: CRS矩阵的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 15:34:57
 */

#ifndef __CHIPSUM_CSR_SERIAL_IMPL_HPP__
#define __CHIPSUM_CSR_SERIAL_IMPL_HPP__

#include <fstream>
#include <vector>
#include <cstdlib>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"

#include "../../common/png_writer.hpp"
#include "../../common/bmp_writer.h"

namespace ChipSum {
namespace Numeric {

template <typename SizeType, typename... Props> struct serial_static_graph {
  ::std::vector<SizeType> row_map;
  ::std::vector<SizeType> col_map;
}__attribute__((aligned));

template <typename ScalarType, typename SizeType, typename... Props>
struct csr_format {
  ::std::vector<ScalarType> vals;
  serial_static_graph<SizeType> graph;
  ::std::size_t col_num;
}__attribute__((aligned));


template <typename ScalarType, typename SizeType, typename... Props>

struct Sparse_Traits<ScalarType, SizeType, SparseTypes::Csr,
                     ChipSum::Backend::Serial, Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Serial,
                             Props...> {

  using sp_type = csr_format<ScalarType, SizeType>;
  using size_type = ::std::size_t;

  using graph_type = serial_static_graph<SizeType>;
  using row_map_type = ::std::vector<SizeType>;
  using col_map_type = ::std::vector<SizeType>;
  using values_type = ::std::vector<ScalarType>;
};

namespace Impl {

namespace Sparse {

template <typename ScalarType, typename SizeType, typename... Props>
// 通过POD数据创建CSR格式矩阵
CHIPSUM_FUNCTION_INLINE void
create(const SizeType nrows, const SizeType ncols, const SizeType annz,
       csr_format<ScalarType, SizeType> &A, SizeType *row_map, SizeType *col_map,
       ScalarType *values) {
  CHIPSUM_UNUSED(ncols);
  A.vals = ::std::vector<ScalarType>(values, values + annz);
  A.graph.row_map = ::std::vector<SizeType>(row_map, row_map + nrows + 1);
  A.graph.col_map = ::std::vector<SizeType>(col_map, col_map + annz);
  A.col_num = ncols;
}

template <typename ScalarType, typename SizeType, typename... Props>
// 创建未初始化的CSR格式矩阵
CHIPSUM_FUNCTION_INLINE void create(csr_format<ScalarType, SizeType> &A,
                                    const SizeType row_map_size,
                                    const SizeType col_map_size) {
  CHIPSUM_UNUSED(row_map_size);
  A.vals.resize(col_map_size);
  A.graph.row_map.resize(col_map_size);
  A.graph.col_map.resize(col_map_size);
}

template <typename ScalarType, typename SizeType, typename... Props>
// b=Ax 一种常用的SpMV接口
CHIPSUM_FUNCTION_INLINE void
mult(::std::size_t M, ::std::size_t N, csr_format<ScalarType, SizeType> &A,
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
// b = beta*b+alpha*A*x 完整的SpMV
CHIPSUM_FUNCTION_INLINE void
mult(ScalarType alpha, csr_format<ScalarType, SizeType> &A,
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
// 命令行打印CSR矩阵的信息。
CHIPSUM_FUNCTION_INLINE void print(csr_format<ScalarType, SizeType> &A,
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
// 命令行打印CSR矩阵的pattern，这个接口适合用来打印些很小的矩阵
CHIPSUM_FUNCTION_INLINE void print_pattern(csr_format<ScalarType, SizeType> &A,
                                   ::std::ostream &out) {
  
  ::std::size_t M = A.graph.row_map.size()-1;
  ::std::size_t N = A.col_num;

  
  ::std::size_t row_entry_cnt = 0;
  ::std::size_t entry_cnt = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = A.graph.row_map[i];
    ::std::size_t end = A.graph.row_map[i+1];
    for (::std::size_t j = 0; j < N; ++j) {
      char info = 'o';
      if(row_entry_cnt < end - start && entry_cnt < A.graph.col_map.size())
      {
      if(A.graph.col_map[start+row_entry_cnt]==j ){
         info = '+';
         ++row_entry_cnt;
         ++entry_cnt;
         
      }
      }
      out<<info<<" ";
    }
    row_entry_cnt = 0;
    out<<::std::endl;
  }
}


template <typename ScalarType, typename SizeType, typename... Props>
// 将稀疏矩阵pattern保存为图片，方便调试和写论文用
CHIPSUM_FUNCTION_INLINE void
save_figure(csr_format<ScalarType,SizeType> &A,
           const char *filename) {

  ::std::size_t M = A.graph.row_map.size() - 1;
  ::std::size_t N = A.col_num;

  ::std::size_t row_entry_cnt = 0;
  ::std::size_t entry_cnt = 0;

  unsigned char *img = static_cast<unsigned char *>(
      ::std::malloc(M * N * 3 * sizeof(unsigned char)));

  char color = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = A.graph.row_map[i];
    ::std::size_t end = A.graph.row_map[i + 1];
    for (::std::size_t j = 0; j < N; ++j) {

      color = 0;
      if (row_entry_cnt < end - start && entry_cnt < A.graph.col_map.size()) {
        if (A.graph.col_map[start + row_entry_cnt] == j) {
          color = 120;
          ++row_entry_cnt;
          ++entry_cnt;
        }
      }

      img[i * N * 3 + j * 3] = color;
      img[i * N * 3 + j * 3 + 1] = color;
      img[i * N * 3 + j * 3 + 2] = color;
    }
    row_entry_cnt = 0;
  }
  ::std::string file_string(filename);

  file_string = &file_string[file_string.find_last_of(".")];

  // ::std::transform(file_string.begin(),
  // file_string.end(),file_string.begin(),::std::toupper); /* 通不过编译，说明我用得不对，暂时还没找到原因. */

  if (file_string == ".bmp" || file_string == ".BMP"/* 补丁写法 */) {
    // 有一些已知的BUG，见用户接口
    ChipSum::Common::flip_bmp(N, M, img);
    ChipSum::Common::write_bmp(N, M, img, filename);
  } else if (file_string == ".png" || file_string == ".PNG"/* 补丁写法 */) {
    ::std::FILE *fp = ::std::fopen(filename, "wb");
    svpng(fp, N, M, img, 0);
    ::std::fclose(fp);
  } else {
    ::std::cerr << "No such format support: " << file_string << endl;
    ::std::cerr << "Saving figure " << filename << " failed!" << endl;
  }
  ::std::free(img);
}


} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
