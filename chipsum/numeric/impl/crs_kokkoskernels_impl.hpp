/*
 * @Description: CRS矩阵的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-10-26 15:36:55
 */

#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__


#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <string>

#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse.hpp>


#include "../../chipsum_macro.h"
#include "../../common/bmp_writer.h"
#include "../../common/png_writer.hpp"
#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"

/*
//  这一部分实现的接口与crs_serial_impl.hpp中类似，不再反复添加注释
*/

static int spm_name = 0;

namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Sparse_Traits<ScalarType, SizeType, SparseTypes::Csr,
                     ChipSum::Backend::KokkosKernels, Props...>
    : public Operator_Traits<ScalarType, SizeType,
                             ChipSum::Backend::KokkosKernels, Props...> {

  using sp_type = KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;
  using size_type = ::std::size_t;

  using graph_type = typename sp_type::staticcrsgraph_type;
  using row_map_type = typename sp_type::row_map_type;
  using col_map_type = typename sp_type::index_type;
  using values_type = typename sp_type::values_type;
};

namespace Impl {

namespace Sparse {

template <typename ScalarType, typename SizeType, typename... Props>
// kokkos实现的创建CSR矩阵
CHIPSUM_FUNCTION_INLINE void
create(const SizeType nrows, const SizeType ncols, const SizeType annz,
       KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
       SizeType *row_map, SizeType *col_map, ScalarType *values) {

  using crs_t = KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  A = crs_t("spm_" + ::std::to_string(spm_name),
            static_cast<typename crs_t::ordinal_type>(nrows),
            static_cast<typename crs_t::ordinal_type>(ncols),
            static_cast<typename crs_t::size_type>(annz), values,
            static_cast<typename crs_t::ordinal_type *>(row_map),
            static_cast<typename crs_t::ordinal_type *>(col_map));
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
create(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
       const ::std::size_t row_map_size, const ::std::size_t col_map_size) {

  using sp_t = KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  typename sp_t::row_map_type row_map(
      "row_map_" + A.label(),
      static_cast<typename sp_t::row_map_type>(row_map_size));
  typename sp_t::entries_type col_map(
      "col_map_" + A.label(),
      static_cast<typename sp_t::entries_type>(col_map_size));

  typename sp_t::staticcrsgraph_type graph(col_map, row_map);

  A = sp_t(A.label(), graph);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
mult(SizeType M, SizeType N,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     const Kokkos::View<ScalarType *> &x, Kokkos::View<ScalarType *> &b) {
  KokkosSparse::spmv("N", static_cast<ScalarType>(1.0), A, x,
                     static_cast<ScalarType>(0.0), b);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
mult(SizeType M, SizeType N, SizeType K,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     const Kokkos::View<ScalarType **> &B, Kokkos::View<ScalarType **> &C) {
  KokkosSparse::spmv("N", static_cast<ScalarType>(1.0), A, B,
                     static_cast<ScalarType>(0.0), C);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
mult(ScalarType alpha,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     Kokkos::View<ScalarType *> &x, ScalarType beta,
     Kokkos::View<ScalarType *> &b) {
  KokkosSparse::spmv("N", alpha, A, x, beta, b);
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
print(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
      ::std::ostream &out) {
  using crs_t =
      typename KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  using row_map_t = typename crs_t::row_map_type::HostMirror;
  using entries_t = typename crs_t::index_type::HostMirror;
  using values_t = typename crs_t::values_type::HostMirror;

  row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
  values_t h_vals = Kokkos::create_mirror_view(A.values);
  entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);

  Kokkos::deep_copy(h_row_map, A.graph.row_map);
  Kokkos::deep_copy(h_vals, A.values);
  Kokkos::deep_copy(h_entries, A.graph.entries);

  out << "spm_" + ::std::to_string(spm_name) << " "
      << "("
      << "rows=" << A.graph.row_map.extent(0) - 1 << ", "
      << "entries=" << h_entries.extent(0) << ")" << ::std::endl;

  out << A.graph.row_map.label() << ": ";
  out << "[";
  for (::std::size_t i = 0; i < h_row_map.extent(0) - 1; ++i) {
    out << h_row_map(i) << ",";
  }
  out << h_row_map(h_row_map.extent(0) - 1) << "]" << ::std::endl;

  out << A.graph.entries.label() << ": ";
  out << "[";
  for (::std::size_t i = 0; i < h_entries.extent(0) - 1; ++i) {
    out << h_entries(i) << ",";
  }
  out << h_entries(h_entries.extent(0) - 1) << "]" << ::std::endl;

  out << A.values.label() << ": ";
  out << "[";
  for (::std::size_t i = 0; i < h_vals.extent(0) - 1; ++i) {
    out << h_vals(i) << ",";
  }
  out << h_vals(h_vals.extent(0) - 1) << "]" << ::std::endl;
}

template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
print_pattern(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
             ::std::ostream &out) {
  using crs_t =
      typename KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  using row_map_t = typename crs_t::row_map_type::HostMirror;
  using entries_t = typename crs_t::index_type::HostMirror;

  row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
  entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);

  Kokkos::deep_copy(h_row_map, A.graph.row_map);
  Kokkos::deep_copy(h_entries, A.graph.entries);

  ::std::size_t M = h_row_map.extent(0) - 1;
  ::std::size_t N = A.numCols();

  ::std::size_t row_entry_cnt = 0;
  ::std::size_t entry_cnt = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = h_row_map[i];
    ::std::size_t end = h_row_map[i + 1];
    for (::std::size_t j = 0; j < N; ++j) {

      char info = 'o';
      if (row_entry_cnt < end - start && entry_cnt < h_entries.extent(0)) {
        if (h_entries[start + row_entry_cnt] == j) {
          info = '+';
          ++row_entry_cnt;
          ++entry_cnt;
        }
      }
      out << info << " ";
    }
    row_entry_cnt = 0;
    out << ::std::endl;
  }
}






template <typename ScalarType, typename SizeType, typename... Props>

CHIPSUM_FUNCTION_INLINE void
save_figure(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
           const char *filename) {
  using crs_t =
      typename KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  using row_map_t = typename crs_t::row_map_type::HostMirror;
  using entries_t = typename crs_t::index_type::HostMirror;

  row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
  entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);

  Kokkos::deep_copy(h_row_map, A.graph.row_map);
  Kokkos::deep_copy(h_entries, A.graph.entries);

  ::std::size_t M = h_row_map.extent(0) - 1;
  ::std::size_t N = A.numCols();

  ::std::size_t row_entry_cnt = 0;
  ::std::size_t entry_cnt = 0;

  unsigned char *img = static_cast<unsigned char *>(
      ::std::malloc(M * N * 3 * sizeof(unsigned char)));


  char color = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = h_row_map[i];
    ::std::size_t end = h_row_map[i + 1];
    for (::std::size_t j = 0; j < N; ++j) {

      color = 0;
      if (row_entry_cnt < end - start && entry_cnt < h_entries.extent(0)) {
        if (h_entries[start + row_entry_cnt] == j) {
          color = 255;
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
    ChipSum::Common::FlipBMP(N, M, img);
    ChipSum::Common::WriteBMP(N, M, img, filename);
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
