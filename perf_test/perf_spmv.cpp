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
#include "../tpls/kokkos-kernels/test_common/KokkosKernels_Test_Structured_Matrix.hpp"
#include "../ChipSum.hpp"


using Scal  = CSFloat;
using Ordinal = default_lno_t;
using Offset  = default_size_type;
using Layout  = default_layout;

int main(int argc, char *argv[]) {

    ChipSum::Common::Init(argc, argv);

    int N = 20;
    if(argc > 1)    N = atoi(argv[1]);

    using device_type = typename Kokkos::Device<
    Kokkos::DefaultExecutionSpace,
    typename Kokkos::DefaultExecutionSpace::memory_space>;
    using mat_type =
    typename KokkosSparse::CrsMatrix<Scal, Ordinal, device_type, void,
    Offset>;
    using values_type = typename mat_type::values_type;

    for(int j=0; j<200; ++j){

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
         mat_structure(0, 0) = N;
         mat_structure(0, 1) = 0;
         mat_structure(0, 2) = 0;
         mat_structure(1, 0) = N;
         mat_structure(1, 1) = 0;
         mat_structure(1, 2) = 0;

         mat_type h_A =
                 Test::generate_structured_matrix2D<mat_type>("FD", mat_structure);

         using row_map_t = typename mat_type::row_map_type::HostMirror;
         using entries_t = typename mat_type::index_type::HostMirror;
         using values_t = typename mat_type::values_type::HostMirror;

         row_map_t h_row_map = Kokkos::create_mirror_view(h_A.graph.row_map);
         values_t h_vals = Kokkos::create_mirror_view(h_A.values);
         entries_t h_entries = Kokkos::create_mirror_view(h_A.graph.entries);

         Kokkos::deep_copy(h_row_map, h_A.graph.row_map);
         Kokkos::deep_copy(h_vals, h_A.values);
         Kokkos::deep_copy(h_entries, h_A.graph.entries);

         auto nrow = h_A.numRows();
         auto ncol = h_A.numCols();
         auto nnz = h_A.nnz();

         CSR A(nrow,ncol,nnz,h_row_map.data(),h_entries.data(),h_vals.data());

         std::vector<Scal> h_v(A.GetRowNum(),1);

         CSVector x(A.GetRowNum(),h_v.data());

         /// \brief 此处用的是A.operator*()，此接口性能不好，但是很方便，
         ///        不用再构造b。如果b已经构造好了，建议调用A.SpMV接口，类
         ///        似接口后续会增添algorithm模块完成。目前先凑合吧。
         CSVector b = A*x;

         int repeat = 100;
         /// \brief 暂时用Kokkos的Timer充数吧
         Kokkos::Timer timer;
         for(int i=0;i<repeat;++i){
             /// \brief 采用此接口免去一个Vector的构造调用，
             ///        性能更优。（6-10x 带宽提升）
             ///        构造函数和析构函数开销真的挺大的。。
             A.SPMV(x,b);
         }
         Kokkos::fence();
         double time = timer.seconds();




         /// \brief 带宽计算公式


         double Gbytes = repeat*1.0e-9*(2*nnz+3*N)/time;

        //  cout<<"---------------------ChipSum Perf Test "<<j+1<<
        //        "---------------------"<<endl;
        // Kokkos::DefaultExecutionSpace::print_configuration(cout,true);
        if(j==0){
            Kokkos::DefaultExecutionSpace::print_configuration(cout,true);
            cout<<"---------------------ChipSum AXPBY SPMV Test "
                "---------------------"<<endl
                <<"nrow/ncol,  nnz,        GFlops :"<<endl;
        }
  
        //  cout<<"CSR matrix info: "<<endl<<
        //        "rows = "<<nrow<<" , columns = "<<ncol;
        //  cout<<" , nnz = "<<nnz<<endl;
        //  cout<<"SpMV performance : "<<Gbytes<<" GFlops"<<endl;
        cout<<setiosflags(ios::left)<<setw(12)<<nrow<<setw(12)<<nnz<<Gbytes<<endl;

         N*=1.05;


    }

    ChipSum::Common::Finalize();
}


