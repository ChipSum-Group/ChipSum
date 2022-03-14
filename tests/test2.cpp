///
/// \file     test2.cpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-12-21
/// \brief    %stuff%
///
#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

typedef Kokkos::vector<double> vector_type;


struct print_functor{
    vector_type _a;
    print_functor(vector_type a)
        :_a(a)
    {

    }

    KOKKOS_INLINE_FUNCTION void operator()(const int i){
        printf("%d",i);
    }
};

#define K 2
#define B 2
#define N 5
#define M 3

int main(int argc,char* argv[])
{
    ChipSum::Common::Init(argc, argv);
    {
        vector_type v(5);
        v(1)=2;

//        Kokkos::parallel_for(5,print_functor(v));
        double *v1 = static_cast<double *>(std::malloc(B*N*N * sizeof(double)));
        for (int i = 0; i < B*N*N ; ++i) {
            v1[i] = double(1);
        }

        double *v2 = static_cast<double *>(std::malloc(B*N*1 * sizeof(double)));
        for (int i = 0; i < B*N*1 ; ++i) {
            v2[i] = double(2);
        }

        double *v3 = static_cast<double *>(std::malloc(B*N*1 * sizeof(double)));
        for (int i = 0; i < B*N*1 ; ++i) {
            v3[i] = double(0);
        }

        Tensor<3> a(B,N,N,v1);
        Tensor<3> a1(B,N,1,v2);
        Tensor<3> a2(B,N,1,v3);
        a.Print();
        // a(0,0,1) = 100;
        std::cout<<a.GetDimthNum(2)<<std::endl;
        a.HostToDevice();
        a.Print();

        a1.Print();
        a2.Print();

        a.GEMV(a1, a2);
        a2.Print();
        
        
        
        double *b1 = static_cast<double *>(std::malloc(K*B*N*M * sizeof(double)));
        for (int i = 0; i < K*B*N*M ; ++i) {
            b1[i] = double(1);
        }
        double *b2 = static_cast<double *>(std::malloc(K*B*M*N * sizeof(double)));
        for (int i = 0; i < K*B*N*M ; ++i) {
            b2[i] = double(2);
        }
        double *b3 = static_cast<double *>(std::malloc(K*B*N*N * sizeof(double)));
        for (int i = 0; i < K*B*N*N ; ++i) {
            b3[i] = double(0);
        }

        Tensor<4> Btensor(K,B,N,M,b1);
        Btensor.Print();
        Tensor<4> Btensor1(K,B,M,N,b2);
        Btensor1.Print();
        Tensor<4> Btensor2(K,B,N,N,b3);
        Btensor2.Print();

        Btensor.GEMM(Btensor1, Btensor2);
        Btensor2.Print();

    }
    ChipSum::Common::Finalize();
}
