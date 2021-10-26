/*
 * @Author: Li Kunyun
 * @Date: 2021-09-07 23:00:38
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-09-24 10:20:58
 * @Description: 
 */
#ifndef __CHIPSUM_DISTRIBUTION_SPARSE_MATRIX_HPP__
#define __CHIPSUM_DISTRIBUTION_SPARSE_MATRIX_HPP__

#include "../numeric/sparse_matrix.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "MPI/mpi_kokkos_impl.hpp"
#endif

namespace ChipSum {
namespace Distribution{
template <typename ScalarType, typename LocalSizeType,typename GlobalSizeType, typename SpFormat,
          typename BackendType, typename... Props>
class Sparse{

    typedef ChipSum::Numeric::SparseMatrix<ScalarType,LocalSizeType,SpFormat,BackendType> NodeSparseType;
    typedef MPIPack<ScalarType,LocalSizeType,GlobalSizeType> MPIPackType;


    NodeSparseType __mat;
    MPIPackType __mpi_pack;
    
    
};
}
}

#endif // __CHIPSUM_DISTRIBUTION_SPARSE_MATRIX_HPP__
