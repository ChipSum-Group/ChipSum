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



struct SparseMatrixBase{};

struct Csr: public SparseMatrixBase{};

}

}
#endif // SPARSEMATRIXTYPES_H
