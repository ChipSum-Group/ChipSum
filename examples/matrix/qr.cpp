/* * * * * * * * * * * * * * * * * * * * *
 *   File:     lu.cpp
 *   Author:   Zhou Xingbin
 *   group:    CDCS-HPC
 *   Time:     2022-02-18
 * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using namespace std;
#include "../ChipSum.hpp"



int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {
        int M = 5;

        CSFloat *m1 = static_cast<CSFloat *>(std::malloc(M*M * sizeof(CSFloat)));

        for(int i=0;i<M*M;++i)
        {
            // A1[i] = CSFloat(rand()) / RAND_MAX;
            m1[i]=i+1;

        }
        CSMatrix A(M, M, m1);
        typedef ChipSum::Numeric::Vector<CSFloat,
        ChipSum::Backend::DefaultBackend> CSVector;
        CSVector a(M);
        CSVector b(M);
        std::cout<<"origin matrix:"<<std::endl;
        A.Print();

        A.QR(a,b);
        std::cout<<"qr tau:"<<std::endl;
        a.Print();
        std::cout<<"qr w:"<<std::endl;
        b.Print();
        std::cout<<"qr matrix:";
        A.Print();


         //CSTensor LU
        int N = 3;
        CSFloat *m2 = static_cast<CSFloat *>(std::malloc(N*M*M*sizeof(CSFloat)));
        for(int i=0;i<N;++i)
            for(int j=0;j<M*M;++j)
                m2[i*M*M+j] = m1[j];
        
        CSTensor<3> A2(N, M, M, m2);
        CSMatrix a2(N, M);
        CSMatrix b2(M, M);
        std::cout<<"origin tensor:"<<std::endl;
        // A2.Print();
        A2.QR(a2,b2);
        std::cout<<"qr tau:"<<std::endl;
        a2.Print();
        std::cout<<"qr w:"<<std::endl;
        b2.Print();
        std::cout<<"qr tensor:"<<std::endl;
        A2.Print();



        
        std::free(m1);
        std::free(m2);
    }
    ChipSum::Common::Finalize();
}