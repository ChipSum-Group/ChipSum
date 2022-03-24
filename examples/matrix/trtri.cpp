/* * * * * * * * * * * * * * * * * * * * *
 *   File:     lu.cpp
 *   Author:   Zhou Xingbin
 *   group:    CDCS-HPC
 *   Time:     2022-03-21
 * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
#include "../ChipSum.hpp"

int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {
        //CSMatrix LU
        int M = 5;
        CSFloat *a = static_cast<CSFloat *>(std::malloc(M*M * sizeof(CSFloat)));
        for(int i=0;i<M;++i)
          for(int j=0;j<M;++j)
        {
            if(j>=i) 
              a[i*M+j] = CSFloat(rand()) / RAND_MAX;
            else
              a[i*M+j] = 0.0;
        }
        CSMatrix A(M, M, a);

        std::cout<<"origin matrix A:"<<std::endl;
        A.Print();
        //最后默认参数"N"
        // A.TRTRI("U","N");
        std::cout<<"return expect 0: "<<A.TRTRI("U")<<std::endl;
        std::cout<<"inverse:A"<<std::endl;
        A.Print();   
        
        std::free(a);

    }
    ChipSum::Common::Finalize();
}