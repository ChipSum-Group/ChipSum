
///
/// \file     densemat_serial_impl.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-05
/// \brief    %stuff%
///
#ifndef __CHIPSUM_DENSEMAT_SERIAL_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_SERIAL_IMPL_HPP__

#include <cassert>
#include <fstream>
#include <vector>

#include "../../../chipsum_macro.h"
#include "../../numeric_traits.hpp"

#include "densemat_serial_gemv_impl.hpp"
#include "densemat_serial_gemm_impl.hpp"
#include "densemat_serial_scal_impl.hpp"


/// 该头文件包含了稠密矩阵的算子实现。
/// 由于一些前期失误，dense matrix被设计成了类似blas的集合接口
/// 似乎是适合于处理小规模稠密矩阵运算的，这样的功能不太具有现实意义
/// 希望后续维护者能进行一些改进，改进建议如下：

/// 我希望dense matrix实际上依旧是面向代数系统求解的，
/// 所以它应当被设计成非常适合LU分解等操作的数据结构。
/// 也就是说矩阵的每一行要抽象为一个向量vector
/// 这样做的目的是为了能兼容如axpby操作，实现可分析扩展的算子组合
/// 例如高斯消元实际上就是不断地做axpby，最终形成三角阵
/// serial实现似乎不需要大规模变动底层数据结构，
/// 但kokkos后端的实现可能需要将view<ValueType**>变为MultiVector，
/// 这部分工作我应该会亲自去实现，但后续的集成工作需要注意兼容。

/* 根据chipsum目前的工作来看，数据结构设计的重要性是远远大于算法的。 */

namespace ChipSum {

namespace Numeric {

template<typename ValueType>
struct serial_densemat{

    ::std::size_t nrow;
    ::std::size_t ncol;
    ::std::vector<ValueType> data;
};

template <typename ValueType,typename... Props>
struct DenseMatrix_Traits<ValueType, ChipSum::Backend::Serial,
        Props...>
        : public Operator_Traits<ValueType> {

    using matrix_type = serial_densemat<ValueType>;

    using value_type =  ValueType;

    using size_type =  ::std::size_t;
};





namespace Impl {
namespace DenseMat {

template <typename ValueType>
// 创建未初始化的矩阵
CHIPSUM_FUNCTION_INLINE void create(
        serial_densemat<ValueType> &A,
        const ::std::size_t M,
        const ::std::size_t N) {

    A.data = ::std::vector<ValueType>(M * N);
    A.nrow = M;
    A.ncol = N;
}

template <typename ValueType>
// 将POD数据填入矩阵
CHIPSUM_FUNCTION_INLINE void create(serial_densemat<ValueType> &A,
                                  const ::std::size_t M,
                                  const ::std::size_t N,
                                  ValueType *src
                                  ) {
    A.data = ::std::vector<ValueType>(src, src + M * N);
    A.nrow = M;
    A.ncol = N;
}





template <typename ValueType>
// 获取A(i,j)
CHIPSUM_FUNCTION_INLINE ValueType &
get_item(serial_densemat<ValueType> &A,
         const std::size_t i,
         const std::size_t j) {

    return A.data[i * A.ncol + j];
}

template <typename ValueType,typename OStreamT>
// 打印矩阵信息
CHIPSUM_FUNCTION_INLINE void print(serial_densemat<ValueType> &A,
                                   OStreamT &out) {

    ::std::size_t M = A.nrow;
    ::std::size_t N = A.ncol;

    cout << "dense_mat_serial("
         << M <<","
         << N <<")"
         << ":" << endl;

    for (::std::size_t i = 0; i < M; ++i) {

        out << " "
            << "[";
        for (::std::size_t j = 0; j < N - 1; ++j) {
            out << A.data[i * N + j] << ", ";
        }
        out << A.data[i * N + N - 1] << "]" << endl;
    }
    out << endl;
}

} // namespace DenseMat
} // namespace Impl
} // namespace Numeric
} // namespace ChipSum

#endif // __CHIPSUM_DENSEMAT_BLAS_IMPL_HPP__
