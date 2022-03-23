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

struct print_functor
{
    vector_type _a;
    print_functor(vector_type a)
        : _a(a)
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(const int i)
    {
        printf("%d", i);
    }
};

#define K 2
#define B 2
#define N 5
#define M 3

int main(int argc, char *argv[])
{
    ChipSum::Common::Init(argc, argv);
    {
        vector_type v(5);
        v(1) = 2;

        //        Kokkos::parallel_for(5,print_functor(v));
        double *v1 = static_cast<double *>(std::malloc(B * N * N * sizeof(double)));
        for (int i = 0; i < B * N * N; ++i)
        {
            v1[i] = double(1);
        }

        double *v2 = static_cast<double *>(std::malloc(B * N * 1 * sizeof(double)));
        for (int i = 0; i < B * N * 1; ++i)
        {
            v2[i] = double(2);
        }

        double *v3 = static_cast<double *>(std::malloc(B * N * 1 * sizeof(double)));
        for (int i = 0; i < B * N * 1; ++i)
        {
            v3[i] = double(0);
        }

        CSTensor<3> tmp(2,2,2);

        CSTensor<3> a(B,N,N,v1);
        CSTensor<3> a1(B,N,1,v2);
        CSTensor<3> a2(B,N,1,v3);
        a.Print();
        a(0,0,1) = 100;
        std::cout<<a.GetDimthNum(2)<<std::endl;
        a.HostToDevice();
        a.Print();

        a1.Print();
        a2.Print();

        a.GEMM(a1, a2);
        a2.Print();

        double *b1 = static_cast<double *>(std::malloc(K * B * N * M * sizeof(double)));
        for (int i = 0; i < K * B * N * M; ++i)
        {
            b1[i] = double(1);
        }
        double *b2 = static_cast<double *>(std::malloc(K * B * M * N * sizeof(double)));
        for (int i = 0; i < K * B * N * M; ++i)
        {
            b2[i] = double(2);
        }
        double *b3 = static_cast<double *>(std::malloc(K * B * N * N * sizeof(double)));
        for (int i = 0; i < K * B * N * N; ++i)
        {
            b3[i] = double(0);
        }

        CSTensor<4> Btensor(K,B,N,M,b1);
        Btensor.Print();
        CSTensor<4> Btensor1(K,B,M,N,b2);
        Btensor1.Print();
        CSTensor<4> Btensor2(K,B,N,N,b3);
        Btensor2.Print();

        // Btensor.GEMM(Btensor1, Btensor2);
        Btensor2.Print();

        double *b4 = static_cast<double *>(std::malloc(K*B*M*1 * sizeof(double)));
        for (int i = 0; i < K*B*M*1 ; ++i) {
            b4[i] = double(4);
        }
        double *b5 = static_cast<double *>(std::malloc(K*B*N*1 * sizeof(double)));
        for (int i = 0; i < K*B*N*1 ; ++i) {
            b5[i] = double(5);
        }
        CSTensor<4> Btensor4(K,B,M,1,b4);
        Btensor4.Print();
        CSTensor<4> Btensor5(K,B,N,1,b5);
        Btensor5.Print();

        Btensor.GEMV(Btensor4, Btensor5);
        Btensor5.Print();

        double *m1 = static_cast<double *>(std::malloc(5*5 * sizeof(double)));
        for (int i = 0; i < 5*5 ; ++i) {
            m1[i] = double(1);
        }
        CSMatrix M1(5,5,m1);

        double *m2 = static_cast<double *>(std::malloc(5 * sizeof(double)));
        for (int i = 0; i < 5 ; ++i) {
            m2[i] = double(i);
        }
        CSVector M2(5,m2);

        M1.SetCol(1,M2);
        M1.SetRow(1,M2);
        M1.Print();

        M1.GetRowCopy(1, M2);
        M2.Print();
        M1.GetColCopy(1, M2);
        M2.Print();
        M1.GetRowSlice(3,0,5,M2); //等大小取，slice大小和取出后存储vector大小相同
        M2.Print();

        CSMatrix M3(2,2);
        M1.GetPartSlice(1,1,3,3,M3);
        M3.Print();
    }
    ChipSum::Common::Finalize();
}
