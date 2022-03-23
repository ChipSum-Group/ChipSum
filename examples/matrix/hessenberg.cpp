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

        CSFloat *A1 = static_cast<CSFloat *>(std::malloc(M*M * sizeof(CSFloat)));

        for(int i=0;i<M*M;++i)
        {
            // A1[i] = CSFloat(rand()) / RAND_MAX;
            A1[i]=i+1;

        }
        CSMatrix A(M, M, A1);
        // typedef ChipSum::Numeric::Vector<CSFloat,
        // ChipSum::Backend::DefaultBackend> CSVector;
        // CSVector a(M);
        // CSVector b(M);
        CSVector a(M);
        CSVector b(M);
        std::cout<<"origin matrix:"<<std::endl;
        A.Print();

        A.HESSENBERG(a,b);
        std::cout<<"tau vector:"<<std::endl;
        a.Print();

        std::cout<<"w:"<<std::endl;
        b.Print();

        std::cout<<"calculate matrix:";
        A.Print();    
        
        std::free(A1);
    }
    ChipSum::Common::Finalize();
}