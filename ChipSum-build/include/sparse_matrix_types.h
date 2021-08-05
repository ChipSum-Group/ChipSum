/* * * * * * * * * * * * * * * * * * * * *
*   File:     sparse_matrix_types.h
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_TYPES_H__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_TYPES_H__


namespace ChipSum {
namespace  Numeric{
namespace  SparseTypes{

struct SparseMatrixBase{};

struct Csr: public SparseMatrixBase{};

struct Csc:public SparseMatrixBase{};

}
}
}
#endif // __CHIPSUM_NUMERIC_SPARSE_MATRIX_TYPES_H__
