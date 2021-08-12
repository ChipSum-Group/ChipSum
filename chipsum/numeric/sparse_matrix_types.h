/*
 * @Description: 稀疏矩阵类型（模板参数）
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:44:32
 */

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
