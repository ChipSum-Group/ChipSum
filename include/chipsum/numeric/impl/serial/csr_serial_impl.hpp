///
/// \file     csr_serial_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-10
/// \brief    %stuff%
///

#ifndef __CHIPSUM_CSR_SERIAL_IMPL_HPP__
#define __CHIPSUM_CSR_SERIAL_IMPL_HPP__

#include <fstream>
#include <vector>
#include <cstdlib>

#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"
#include "../../sparse_matrix_types.h"

#include "../../../common/png_writer.hpp"
#include "../../../common/bmp_writer.h"


#include "csr_serial_spmv_impl.hpp"


namespace ChipSum {
namespace Numeric {


struct serial_static_graph {
  ::std::vector<::std::size_t> row_map;
  ::std::vector<::std::size_t> col_map;
}__attribute__((aligned));

template <typename ValueType>
struct serial_csr_format {
  ::std::vector<ValueType> vals;
  serial_static_graph graph;
  ::std::size_t col_num;
}__attribute__((aligned));


template <typename ValueType,typename ...Props>

struct Sparse_Traits<ValueType,
        ChipSum::Backend::Serial,
        ChipSum::Numeric::SparseTypes::Csr,
        Props...>
    : public Operator_Traits<ValueType> {

  using sp_type = serial_csr_format<ValueType>;
  using size_type = ::std::size_t;
    using ordinal_type = ::std::size_t;
  using backend_type = ChipSum::Backend::Serial;
    using format_type = ChipSum::Numeric::SparseTypes::Csr;

  using graph_type = serial_static_graph;
  using row_map_type = ::std::vector<::std::size_t>;
  using col_map_type = ::std::vector<::std::size_t>;
  using values_type = ::std::vector<ValueType>;
    using value_type = ValueType;
};

namespace Impl {

namespace Sparse {

template <typename ValueType,typename S1,typename S2,typename S3,typename S4,typename S5 >
// 通过POD数据创建CSR格式矩阵
CHIPSUM_FUNCTION_INLINE void
create(serial_csr_format<ValueType> &A,
       S1 nrows,
       S2 ncols,
       S3 annz,
       S4 *row_map,
       S5 *col_map,
       ValueType *values) {
//  CHIPSUM_UNUSED(ncols);



  A.vals = ::std::vector<ValueType>(values, values + annz);
  A.graph.row_map = ::std::vector<::std::size_t>(row_map, row_map + nrows + 1);
  A.graph.col_map = ::std::vector<::std::size_t>(col_map, col_map + annz);
  A.col_num = ncols;
}

template <typename ValueType>
// 创建未初始化的CSR格式矩阵
CHIPSUM_FUNCTION_INLINE void create(serial_csr_format<ValueType> &A,
                                    const ::std::size_t nrows,
                                    const ::std::size_t ncols,
                                    const ::std::size_t annz) {

  A.vals.resize(annz);
  A.graph.row_map.resize(nrows+1);
  A.graph.col_map.resize(annz);
  A.col_num = ncols;
}



template <typename ValueType,typename OStreamT>
// 命令行打印CSR矩阵的信息。
CHIPSUM_FUNCTION_INLINE void print(serial_csr_format<ValueType> &A,
                                   OStreamT &out) {
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

template <typename ValueType,typename OStreamT>
// 命令行打印CSR矩阵的pattern，这个接口适合用来打印些很小的矩阵
CHIPSUM_FUNCTION_INLINE void print_pattern(serial_csr_format<ValueType> &A,
                                   OStreamT &out) {
  
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


template <typename ValueType>
// 将稀疏矩阵pattern保存为图片，方便调试和写论文用
CHIPSUM_FUNCTION_INLINE void
save_figure(serial_csr_format<ValueType> &A,
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
          color = static_cast<char>(250);
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
