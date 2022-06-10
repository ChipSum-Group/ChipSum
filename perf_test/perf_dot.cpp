/* * * * * * * * * * * * * * * * * * * * *
 *   File:     test.cpp
 *   Author:   Li Kunyun
 *   group:    CDCS-HPC
 *   Time:     2021-07-28
 * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include "../ChipSum.hpp"
#include <KokkosKernels_IOUtils.hpp>



int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {   int N=100;

        // cnt记录flops两次差值比例在1%以内次数，flag记录连续性，连续5次为true
        // 连续5次差值比例在1%以内，break
        int cnt = 0;
        bool flag = false;
        double pre = 0;

        for(int j=0; j<200; ++j){
            CSFloat *v1 = static_cast<CSFloat *>(std::malloc(N * sizeof(CSFloat)));
            CSFloat *v2 = static_cast<CSFloat *>(std::malloc(N * sizeof(CSFloat)));
            
            for (int i = 0; i < N; ++i) {
                v1[i] = CSFloat(i);
                v2[i] = CSFloat(i);
            }
            
            CSVector a(N,v1); // a = {0,1,2,3,4...}
            CSVector b(N,v2); // b = {0,1,2,3,4...}

            // Scalar r ;
            int repeat = 100;
            
            Kokkos::Timer timer;
            for(int i=0;i<repeat;++i){
            a.Dot(b);
            }
            Kokkos::fence();
            double time = timer.seconds();
            

            /// \brief 带宽计算公式
            double Gbytes = repeat*1.0e-9*(2.0*N-1)/time;

            cout << "CSVector size: " << setiosflags(ios::left)<<setw(12)<<N<<"GFlops :"<<Gbytes<<endl;
             

            N*=1.1;          
            std::free(v1);
            std::free(v2);

            if(abs(Gbytes-pre)/pre < 0.01){
                cnt += 1;
                flag = true;
                // cout<<"cnt: "<< cnt<<endl;
            }
            else{
                cnt = 0;
                flag = false;
            }
            pre = Gbytes;

            if(cnt==5 && flag) break;  
        }
    }
    ChipSum::Common::Finalize();
}
