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
        int M = 10;

        CSFloat *A1 = static_cast<CSFloat *>(std::malloc(M*M * sizeof(CSFloat)));

        for(int i=0;i<M*M;++i)
        {
            A1[i] = CSFloat(rand()) / RAND_MAX;
        }
        Matrix A(M, M, A1);
        std::cout<<"origin matrix:"<<std::endl;
        A.Print();

        A.LU();
        std::cout<<"lu matrix:"<<std::endl;
        A.Print();      
        
        std::free(A1);
    }
    ChipSum::Common::Finalize();
}