/*
 * @Description: CRS矩阵的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:45:41
 */


#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__

#include <KokkosKernels_Handle.hpp>
#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosSparse_spmv.hpp>

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"

static int spm_name = 0;

namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Sparse_Traits<ScalarType, SizeType, SparseTypes::Csr,
                     ChipSum::Backend::KokkosKernels, Props...>
    : public Operator_Traits<ScalarType, SizeType,
                             ChipSum::Backend::KokkosKernels, Props...> {

  using sp_type = KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;
  using size_type = std::size_t;

  using graph_type = typename sp_type::staticcrsgraph_type;
  using row_map_type = typename sp_type::row_map_type;
  using col_map_type = typename sp_type::index_type;
  using values_type = typename sp_type::values_type;
};

namespace Impl {

namespace Sparse {

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Fill：利用POD数据对稀疏矩阵进行填充
 * @param A：稀疏矩阵
 * @param nrows：行数
 * @param ncols：列数
 * @param annz：非零元数
 * @param row_map：行向量
 * @param col_map：entry向量
 * @param values：非零元向量
 */
CHIPSUM_FUNCTION_INLINE void
Create(const SizeType nrows, const SizeType ncols, const SizeType annz,
       KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
       SizeType *row_map, SizeType *col_map, ScalarType *values) {

  using crs_t = KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  A = crs_t("spm_" + std::to_string(spm_name),
            static_cast<typename crs_t::ordinal_type>(nrows),
            static_cast<typename crs_t::ordinal_type>(ncols),
            static_cast<typename crs_t::size_type>(annz), values,
            static_cast<typename crs_t::ordinal_type *>(row_map),
            static_cast<typename crs_t::ordinal_type *>(col_map));
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Create
 * @param A
 * @param row_map_size
 * @param col_map_size
 */
CHIPSUM_FUNCTION_INLINE void
Create(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
       const std::size_t row_map_size, const std::size_t col_map_size) {

  using sp_t = KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device>;

  typename sp_t::row_map_type row_map(
      "row_map_" + A.label(),
      static_cast<typename sp_t::row_map_type>(row_map_size));
  typename sp_t::entries_type col_map(
      "col_map_" + A.label(),
      static_cast<typename sp_t::entries_type>(col_map_size));

  typename sp_t::staticcrsgraph_type graph(col_map, row_map);

  A = spt(A.label(), graph);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Mult b=Ax
 * @param A
 * @param x
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void
Mult(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     const Kokkos::View<ScalarType *> &x, Kokkos::View<ScalarType *> &b) {
  KokkosSparse::spmv("N", static_cast<ScalarType>(1.0), A, x,
                     static_cast<ScalarType>(0.0), b);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Mult
 * @param A
 * @param x
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void
Mult(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     const Kokkos::View<ScalarType **> &B, Kokkos::View<ScalarType **> &C) {
  KokkosSparse::spmv("N", static_cast<ScalarType>(1.0), A, B,
                     static_cast<ScalarType>(0.0), C);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @brief Mult
 * @param alpha
 * @param A
 * @param x
 * @param beta
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void
Mult(ScalarType alpha,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     Kokkos::View<ScalarType *> &x, ScalarType beta,
     Kokkos::View<ScalarType *> &b) {
  KokkosSparse::spmv("N", alpha, A, x, beta, b);
}

} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
