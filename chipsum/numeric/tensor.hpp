///
/// \file     tensor.hpp
/// \author   Yukai Xiang
/// \group    CDCS-HPC
/// \date     2022-03-07
/// \brief    张量用户接口
///

#ifndef __CHIPSUM_TENSOR_HPP__
#define __CHIPSUM_TENSOR_HPP__


#include "impl/kokkoskernels/tensor_kokkoskernels_impl.hpp"

#include "impl/serial/densemat_serial_impl.hpp"

#include "numeric_traits.hpp"
#include "scalar.hpp"
#include "vector.hpp"
#include "dense_matrix.hpp"

namespace ChipSum {
namespace Numeric {


template <size_t NDIM, typename... Props>
class TensorNDIM {

public:
    using traits = Tensor_Traits<NDIM, Props...>;
    
    using tensor_type = typename traits::tensor_type;
    using tensor_type_ref = typename std::add_lvalue_reference<tensor_type>::type;
    using const_tensor_type_ref = typename std::add_const<tensor_type_ref>::type;

    using size_type = typename traits::size_type;
    using const_size_type = const size_type;
    using const_size_type_ref = const size_type&;

    using value_type = typename traits::value_type;
    using backend_type = typename traits::backend_type;

    using vector_type = ChipSum::Numeric::Vector< Props...>;
    using matrix_type =  ChipSum::Numeric::DenseMatrix< Props...>;

private:
    tensor_type __data;
    size_type __nnum;
    size_type __nbatch;
    size_type __nrow;
    size_type __ncol;

public:

    ///
    /// \brief TensorNDIM 构造一个3维张量，该张量未初始化
    /// \param Batch Batch数
    /// \param M M 行数
    /// \param N N 列数
    ///
    CHIPSUM_DECLARED_FUNCTION TensorNDIM(const_size_type Batch, const_size_type M, const_size_type N)
        : __nbatch(Batch), __nrow(M), __ncol(N) {
        ChipSum::Numeric::Impl::Tensor::create(__data, Batch, M, N);
    }

    ///
    /// \brief TensorNDIM 构造一个3维张量，并张量初始化
    /// \param Batch Batch数
    /// \param M M 行数
    /// \param N N 列数
    ///
    CHIPSUM_DECLARED_FUNCTION TensorNDIM(const_size_type Batch, const_size_type M, const_size_type N, value_type *src)
        : __nbatch(Batch), __nrow(M), __ncol(N) {
        // ChipSum::Numeric::Impl::Tensor::create(__data, Batch, M, N, src);
        ChipSum::Numeric::Impl::Tensor::create(src, __data, Batch, M, N);
    }

    ///
    /// \brief TensorNDIM 构造一个4维张量，该张量未初始化
    /// \param Num 个数
    /// \param Batch Batch数
    /// \param M M 行数
    /// \param N N 列数
    ///
    CHIPSUM_DECLARED_FUNCTION TensorNDIM(const_size_type Num, const_size_type Batch, const_size_type M, const_size_type N)
        : __nnum(Num), __nbatch(Batch), __nrow(M), __ncol(N) {
        ChipSum::Numeric::Impl::Tensor::create(__data, Num, Batch, M, N);
    }

    ///
    /// \brief TensorNDIM 构造一个4维张量，并张量初始化
    /// \param Num 个数
    /// \param Batch Batch数
    /// \param M M 行数
    /// \param N N 列数
    ///
    CHIPSUM_DECLARED_FUNCTION TensorNDIM(const_size_type Num, const_size_type Batch, const_size_type M, const_size_type N, value_type *src)
        : __nnum(Num), __nbatch(Batch), __nrow(M), __ncol(N) {
        ChipSum::Numeric::Impl::Tensor::create(src, __data, Num, Batch, M, N);
    }


    ///
    /// \brief GetData 获取Tensor数据
    /// \return 后端数据
    ///
    CHIPSUM_FUNCTION_INLINE tensor_type_ref GetData() {
        return __data;
    }

    ///
    /// \brief GetDimthNum 获取指定维度的大小
    /// \return 后端数据
    ///
    CHIPSUM_FUNCTION_INLINE size_type GetDimthNum(size_t Dimth) { return __data.extent(Dimth); }

    ///
    /// \brief operator () 获取host端A(i,j,k)
    /// \param i 第0维度索引
    /// \param j 第1维度索引
    /// \param k 第2维度索引
    /// \return A(i,j,k)
    ///
    template<typename ...Args>
    CHIPSUM_FUNCTION_INLINE value_type &operator()(Args ...args) {
        return __data.h_view(args...);
    }

    ///
    /// \brief operator () 获取host端A(i,j,k)（只读）
    /// \param i 第0维度索引
    /// \param j 第1维度索引
    /// \param k 第2维度索引
    /// \return A(i,j,k)
    ///
    template<typename ...Args>
    CHIPSUM_FUNCTION_INLINE const value_type &operator()(Args ...args) const{
        return __data.h_view(args...);
    }
    

    ///
    /// \brief At 获取device端A(i,j,k)
    /// \param i 第0维度索引
    /// \param j 第1维度索引
    /// \param k 第2维度索引
    /// \return A(i,j,k)
    ///
    template<typename ...Args>
    CHIPSUM_FUNCTION_INLINE value_type At(Args ...args) {
        return __data.d_view(args...);
    }

    ///
    /// \brief Device端到Host端数据深拷贝
    ///
    CHIPSUM_FUNCTION_INLINE void DeviceToHost(){
        ChipSum::Numeric::Impl::Tensor::device_to_host(__data);
    }

    ///
    /// \brief Host端到Device端数据深拷贝
    ///
    CHIPSUM_FUNCTION_INLINE void HostToDevice(){
        ChipSum::Numeric::Impl::Tensor::host_to_device(__data);
    }


    ///
    /// \brief GEMM C=A*B 当C为已初始化的矩阵时，强烈建议采用此接口进行GEMM运算
    /// \param B  参与运算的另一矩阵
    /// \param C  结果
    ///
    CHIPSUM_FUNCTION_INLINE void GEMM(TensorNDIM &m, TensorNDIM &ret) {
        ChipSum::Numeric::Impl::Tensor::batch_gemm(__data, m.GetData(), ret.GetData());
    }

    ///
    /// \brief GEMV y=A*x 
    /// \param x  参与运算的向量
    /// \param y  结果
    ///
    CHIPSUM_FUNCTION_INLINE void GEMV(TensorNDIM &x, TensorNDIM& y) {
        ChipSum::Numeric::Impl::Tensor::batch_gemv(__data, x.GetData(), y.GetData());
    }

    ///
    /// \brief LU分解，结果存入原矩阵
    /// \param tiny 分解精度,默认为0
    ///
    CHIPSUM_FUNCTION_INLINE void LU(const value_type tiny = 0) {
        ChipSum::Numeric::Impl::Tensor::batch_lu(__data,tiny);
    }

    ///
    /// \brief QR分解，结果存入原矩阵
    /// \param 输出矩阵系数tau
    /// \param 输出矩阵系数w
    ///
    CHIPSUM_FUNCTION_INLINE void QR(matrix_type &x, matrix_type& y) {
        ChipSum::Numeric::Impl::Tensor::batch_qr(__data,x.GetData(),y.GetData());
    }

    template<typename OStreamT=std::ostream>
    ///
    /// \brief Print 打印函数
    /// \param out 输出流
    ///
    CHIPSUM_FUNCTION_INLINE void Print(OStreamT &out = std::cout) {
        ChipSum::Numeric::Impl::Tensor::print(
                    __data,out);
    }

};

} // End namespace Numeric
} // End namespace ChipSum

///
/// \brief Tensor 默认的
///
template<size_t NDIM> using CSTensor = ChipSum::Numeric::TensorNDIM<NDIM, CSFloat, ChipSum::Backend::DefaultBackend>;

#endif // __CHIPSUM_TENSOR_HPP__
