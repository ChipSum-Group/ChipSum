/*
 * @Description: CRS矩阵的KokkosKernels实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-17 15:33:26
 */

#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__



#include <KokkosSparse.hpp>
#include <KokkosKernels_default_types.hpp>
#include <fstream>


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
 * @description: 创建未初始化的CRS矩阵
 * @param {*} A 稀疏矩阵（out）
 * @param {*} row_map_size 行邻接表长度
 * @param {*} col_map_size 列邻接表长度
 * @return {*} 
 * @author: Li Kunyun
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

  A = sp_t(A.label(), graph);
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
Mult(SizeType M, SizeType N,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     const Kokkos::View<ScalarType *> &x, Kokkos::View<ScalarType *> &b) {
  KokkosSparse::spmv("N", static_cast<ScalarType>(1.0), A, x,
                     static_cast<ScalarType>(0.0), b);
}

template <typename ScalarType, typename SizeType, typename... Props>
/**
 * @description: SpGEMM（稀疏X稠密）
 * @param {*} M A/C的行数
 * @param {*} N B/C的列数
 * @param {*} K A的列数，B的行数
 * @param {*} A 稀疏矩阵
 * @param {*} B 稠密矩阵
 * @param {*} C 稠密矩阵（out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Mult(SizeType M, SizeType N, SizeType K,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     const Kokkos::View<ScalarType **> &B, Kokkos::View<ScalarType **> &C) {
  KokkosSparse::spmv("N", static_cast<ScalarType>(1.0), A, B,
                     static_cast<ScalarType>(0.0), C);
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
Mult(ScalarType alpha,
     KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
     Kokkos::View<ScalarType *> &x, ScalarType beta,
     Kokkos::View<ScalarType *> &b) {
  KokkosSparse::spmv("N", alpha, A, x, beta, b);
}



template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @description: 打印出稀疏矩阵的数据信息，多用于调试
 * @param {*} A 稀疏矩阵
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
Print(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
      std::ostream &out) 
{
  using crs_t = typename KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>;
  
  using row_map_t = typename crs_t::row_map_type::HostMirror;
  using entries_t = typename crs_t::index_type::HostMirror;
  using values_t = typename crs_t::values_type::HostMirror;

  row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
  values_t h_vals = Kokkos::create_mirror_view(A.values);
  entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);



  Kokkos::deep_copy(h_row_map,A.graph.row_map);
  Kokkos::deep_copy(h_vals,A.values);
  Kokkos::deep_copy(h_entries,A.graph.entries);

  
  
  out<<"spm_"+std::to_string(spm_name)<<" "
  <<"("<<"rows="<<A.graph.row_map.extent(0)-1<<", "
  <<"entries="<<h_entries.extent(0)<<")"<<std::endl;



   out<<A.graph.row_map.label()<<": ";
  out<<"[";
  for(std::size_t i=0;i<h_row_map.extent(0)-1;++i){
     out<<h_row_map(i)<<",";
  }
  out<<h_row_map(h_row_map.extent(0)-1)<<"]"<<std::endl;


  
  out<<A.graph.entries.label()<<": ";
  out<<"[";
  for(std::size_t i=0;i<h_entries.extent(0)-1;++i){
     out<<h_entries(i)<<",";
  }
  out<<h_entries(h_entries.extent(0)-1)<<"]"<<std::endl;
  

  out<<A.values.label()<<": ";
  out<<"[";
  for(std::size_t i=0;i<h_vals.extent(0)-1;++i){
     out<<h_vals(i)<<",";
  }
  out<<h_vals(h_vals.extent(0)-1)<<"]"<<std::endl;
}


template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @description: 打印出稀疏矩阵的数据信息，多用于调试
 * @param {*} A 稀疏矩阵
 * @param {*} out 输出流（in/out）
 * @return {*}
 * @author: Li Kunyun
 */
CHIPSUM_FUNCTION_INLINE void
PrintPattern(KokkosSparse::CrsMatrix<ScalarType, SizeType, default_device> &A,
      std::ostream &out) 
{
  using crs_t = typename KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>;
  
  using row_map_t = typename crs_t::row_map_type::HostMirror;
  using entries_t = typename crs_t::index_type::HostMirror;

  row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
  entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);

  Kokkos::deep_copy(h_row_map,A.graph.row_map);
  Kokkos::deep_copy(h_entries,A.graph.entries);

  std::size_t M = h_row_map.extent(0)-1;
  std::size_t N = A.numCols();


  std::size_t entry_cnt = 0;
   
  for (std::size_t i = 0; i < M; ++i) {
    std::size_t start = h_row_map[i];
    for (std::size_t j = 0; j < N; ++j) {
      char info = 'o';
      if(h_entries[start+entry_cnt]==j){
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
