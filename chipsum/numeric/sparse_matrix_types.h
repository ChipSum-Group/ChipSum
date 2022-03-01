///
/// \file     sparse_matrix_types.h
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-10-27
///
#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_TYPES_H__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_TYPES_H__



namespace ChipSum {
namespace  Numeric{
namespace  SparseTypes{

struct SparseMatrixBase{};

struct Csr: public SparseMatrixBase{};

struct Csc:public SparseMatrixBase{};

struct Coo: public SparseMatrixBase{};

}
}
}
#endif // __CHIPSUM_NUMERIC_SPARSE_MATRIX_TYPES_H__
