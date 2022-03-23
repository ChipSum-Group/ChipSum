/* * * * * * * * * * * * * * * * * * * * *
 *   File:     test.cpp
 *   Author:   Li Kunyun
 *   group:    CDCS-HPC
 *   Time:     2021-07-28
 * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using namespace std;

#include <vector>
#include "../ChipSum.hpp"



int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {   int N=100;
        CSFloat *v1 = static_cast<CSFloat *>(std::malloc(N * sizeof(CSFloat)));
        CSFloat *v2 = static_cast<CSFloat *>(std::malloc(N * sizeof(CSFloat)));
        
        for (int i = 0; i < N; ++i) {
            v1[i] = CSFloat(i);
            v2[i] = CSFloat(i);
        }
        
        CSVector a(N,v1); // a = {0,1,2,3,4...}
        CSVector b(N,v2); // b = {0,1,2,3,4...}

        a.Dot(b);
        Kokkos::fence();
                
        std::free(v1);
        std::free(v2);
    }
    ChipSum::Common::Finalize();
}
