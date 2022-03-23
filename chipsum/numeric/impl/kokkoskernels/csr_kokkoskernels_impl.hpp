///
/// \file     csr_kokkoskernels_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-05
/// \brief    %stuff%
///

#ifndef __CHIPSUM_CSR_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CSR_KOKKOSKERNELS_IMPL_HPP__


#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <string>

#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse.hpp>


#include "../../../chipsum_macro.h"
#include "../../../common/bmp_writer.h"
#include "../../../common/png_writer.hpp"

#include "../../numeric_traits.hpp"
#include "../../sparse_matrix_types.h"

#include "csr_kokkoskernels_spmv_impl.hpp"
#include "csr_kokkoskernels_spgemm_impl.hpp"
#include "csr_kokkoskernels_spilu_impl.hpp"
#include "csr_kokkoskernels_sptrsv_impl.hpp"

/*
//  这一部分实现的接口与crs_serial_impl.hpp中类似，不再反复添加注释
*/

static int spm_name = 0;

namespace ChipSum {
namespace Numeric {

template <typename ValueType,typename OrdinalType,typename SizeType>
struct Sparse_Traits<ValueType,
        ChipSum::Backend::KokkosKernels,
        ChipSum::Numeric::SparseTypes::Csr,
        OrdinalType,
        SizeType>
        : public Operator_Traits<ValueType> {

    using sp_type = KokkosSparse::CrsMatrix<ValueType, OrdinalType, default_device,void,SizeType>;

    using ordinal_type = typename sp_type::ordinal_type;
    using size_type = typename sp_type::size_type;
    using value_type = typename sp_type::value_type;


    using backend_type = ChipSum::Backend::KokkosKernels;
    using format_type = ChipSum::Numeric::SparseTypes::Csr;


    using graph_type = typename sp_type::staticcrsgraph_type;
    using row_map_type = typename sp_type::row_map_type;
    using col_map_type = typename sp_type::index_type;
    using values_type = typename sp_type::values_type;
};

namespace Impl {

namespace Sparse {


#define __matrix_type KokkosSparse::CrsMatrix \
    <ValueType, OrdinalType, default_device,void,SizeType>

template <typename ValueType,
          typename OrdinalType,
          typename SizeType,
          typename S1,
          typename S2,
          typename S3,
          typename S4,
          typename S5 >
// kokkos实现的创建CSR矩阵
CHIPSUM_FUNCTION_INLINE void
create(__matrix_type &A,
       const S1& nrows,
       const S2& ncols,
       const S3& annz,
       S4 *row_map,
       S5 *col_map,
       ValueType *values) {

    A = __matrix_type("spm_" + ::std::to_string(spm_name++),
                    static_cast<typename __matrix_type::ordinal_type>(nrows),
                    static_cast<typename __matrix_type::ordinal_type>(ncols),
                    static_cast<typename __matrix_type::size_type>(annz), values,
                    static_cast<typename __matrix_type::size_type *>(row_map),
                    static_cast<typename __matrix_type::ordinal_type *>(col_map));
}

template <typename ValueType,
          typename OrdinalType,
          typename SizeType>

CHIPSUM_FUNCTION_INLINE void
create(__matrix_type &A,
       const OrdinalType row_map_size, const OrdinalType col_map_size, const OrdinalType) {

    typename __matrix_type::row_map_type row_map(
                "row_map_" + A.label(),
                static_cast<typename __matrix_type::row_map_type>(row_map_size));
    typename __matrix_type::entries_type col_map(
                "col_map_" + A.label(),
                static_cast<typename __matrix_type::entries_type>(col_map_size));

    typename __matrix_type::staticcrsgraph_type graph(col_map, row_map);

    A = __matrix_type(A.label(), graph);
}

template <typename ValueType,
          typename OrdinalType,
          typename SizeType>

CHIPSUM_FUNCTION_INLINE void
print(__matrix_type &A,
      ::std::ostream &out) {
    using crs_t =
    typename __matrix_type;

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

template <typename ValueType,
          typename OrdinalType,
          typename SizeType>

CHIPSUM_FUNCTION_INLINE void
print_pattern(__matrix_type &A,
              ::std::ostream &out) {
    using crs_t =
    typename __matrix_type;

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






template <typename ValueType,
          typename OrdinalType,
          typename SizeType>

CHIPSUM_FUNCTION_INLINE void
save_figure(__matrix_type &A,
            const char *filename) {
    using crs_t =
    typename __matrix_type;

    using row_map_t = typename crs_t::row_map_type::HostMirror;
    using entries_t = typename crs_t::index_type::HostMirror;

    row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
    entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);

    Kokkos::deep_copy(h_row_map, A.graph.row_map);
    Kokkos::deep_copy(h_entries, A.graph.entries);

    SizeType M = h_row_map.extent(0) - 1;
    SizeType N = A.numCols();

    SizeType row_entry_cnt = 0;
    decltype (h_entries.extent(0)) entry_cnt = 0;

    unsigned char *img = static_cast<unsigned char *>(
                ::std::malloc(M * N * 3 * sizeof(unsigned char)));


    char color = 0;

    for (SizeType i = 0; i < M; ++i) {
        SizeType start = h_row_map[i];
        SizeType end = h_row_map[i + 1];
        for (SizeType j = 0; j < N; ++j) {

            color = 0;
            if (row_entry_cnt < end - start && entry_cnt < h_entries.extent(0)) {
                if (h_entries[start + row_entry_cnt] == j) {
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
        ::std::cerr << "No such format support: " << file_string << ::std::endl;
        ::std::cerr << "Saving figure " << filename << " failed!" << ::std::endl;
    }
    ::std::free(img);
}




} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum
#endif

#endif // __CHIPSUM_CSR_KOKKOSKERNELS_IMPL_HPP__
