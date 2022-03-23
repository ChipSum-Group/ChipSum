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
        //CSMatrix LU
        int M = 5;
        CSFloat *a = static_cast<CSFloat *>(std::malloc(M*M * sizeof(CSFloat)));
        for(int i=0;i<M*M;++i)
        {
            a[i] = CSFloat(rand()) / RAND_MAX;
        }
        CSMatrix A(M, M, a);
        std::cout<<"origin matrix:"<<std::endl;
        A.Print();
        A.LU();
        std::cout<<"lu matrix:"<<std::endl;
        A.Print();      
        

        //CSTensor LU
        int N = 3;
        CSFloat *a2 = static_cast<CSFloat *>(std::malloc(N*M*M*sizeof(CSFloat)));
        for(int i=0;i<N;++i)
            for(int j=0;j<M*M;++j)
                a2[i*M*M+j] = a[j];
        
        CSTensor<3> A2(N,M,M,a2);
        std::cout<<"origin tensor:"<<std::endl;
        A2.Print();
        A2.LU(0.0);
        std::cout<<"lu tensor:"<<std::endl;
        A2.Print();  

        std::free(a);
        std::free(a2);

    }
    ChipSum::Common::Finalize();
}