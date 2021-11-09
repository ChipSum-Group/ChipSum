///
/// \file     perf_test.cpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-08
/// \brief    %stuff%
///

#include <iostream>
using namespace std;

#include <KokkosKernels_default_types.hpp>

#include "tpl/kokkos-kernels/test_common/KokkosKernels_Test_Structured_Matrix.hpp"
#include "ChipSum.hpp"


using Scal  = default_scalar;
using Ordinal = default_lno_t;
using Offset  = default_size_type;
using Layout  = default_layout;

int main(int argc, char *argv[]) {

    ChipSum::Common::Init(argc, argv);

    using device_type = typename Kokkos::Device<
    Kokkos::DefaultExecutionSpace,
    typename Kokkos::DefaultExecutionSpace::memory_space>;
    using mat_type =
    typename KokkosSparse::CrsMatrix<Scal, Ordinal, device_type, void,
    Offset>;
    using values_type = typename mat_type::values_type;

    {

        /// \brief 创造了一个经典有限差分问题的稀疏矩阵
        ///        大概长这样：
        ///        +  +  o  +  o  o  o  o  o
        ///        +  +  +  o  +  o  o  o  o
        ///        o  +  +  o  o  +  o  o  o
        ///        +  o  o  +  +  o  +  o  o
        ///        o  +  o  +  +  +  o  +  o
        ///        o  o  +  o  +  +  o  o  +
        ///        o  o  o  +  o  o  +  +  o
        ///        o  o  o  o  +  o  +  +  +
        ///        o  o  o  o  o  +  o  +  +
        Kokkos::View<Ordinal* [3], Kokkos::HostSpace> mat_structure(
                    "Matrix Structure", 2);
        mat_structure(0, 0) = 50;
        mat_structure(0, 1) = 0;
        mat_structure(0, 2) = 0;
        mat_structure(1, 0) = 50;
        mat_structure(1, 1) = 0;
        mat_structure(1, 2) = 0;

        mat_type A =
                Test::generate_structured_matrix2D<mat_type>("FD", mat_structure);

        using row_map_t = typename mat_type::row_map_type::HostMirror;
        using entries_t = typename mat_type::index_type::HostMirror;
        using values_t = typename mat_type::values_type::HostMirror;

        row_map_t h_row_map = Kokkos::create_mirror_view(A.graph.row_map);
        values_t h_vals = Kokkos::create_mirror_view(A.values);
        entries_t h_entries = Kokkos::create_mirror_view(A.graph.entries);

        Kokkos::deep_copy(h_row_map, A.graph.row_map);
        Kokkos::deep_copy(h_vals, A.values);
        Kokkos::deep_copy(h_entries, A.graph.entries);

        auto nrow = A.numRows();
        auto ncol = A.numCols();
        auto nnz = A.nnz();

        CSR spm(nrow,ncol,nnz,h_row_map.data(),h_entries.data(),h_vals.data());

        std::vector<Scal> h_v(spm.GetRowNum(),1);

        Vector x(h_v.data(),spm.GetRowNum());

        /// \brief 此处用的是spm.operator*()，此接口性能不好，但是很方便，
        ///        不用再构造b。如果b已经构造好了，建议调用spm.SpMV接口，类
        ///        似接口后续会增添algorithm模块完成。目前先凑合吧。
        Vector b = spm*x;

        int repeat = 100;
        /// \brief 暂时用Kokkos的Timer充数吧
        Kokkos::Timer timer;
        for(int i=0;i<repeat;++i){
            /// \brief 采用此接口免去一个Vector的构造调用，
            ///        性能更优。（6-10x spmv带宽提升）
            ///        构造函数和析构函数开销真的挺大的。。
            spm.SpMV(x,b);
        }
        double time = timer.seconds();




        /// \brief 带宽计算公式
        double Gbytes = repeat*1.0e-9*(sizeof(Ordinal)*nrow +
                                sizeof(Ordinal)*nnz +
                                sizeof(Scal)*nnz +
                                sizeof(Scal)*ncol*2)/time;

        /// \brief Pascal61 (MX350)有46——50 GFlops左右的吞吐量
        ///        这个数据应算不错了。
        cout<<Gbytes<<" GFlops"<<endl;
    }

    ChipSum::Common::Finalize();
}


