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
#include <cstdlib>
#include <ctime>
#include <KokkosKernels_IOUtils.hpp>



int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {
        int M = 10;

        // cnt记录flops两次差值比例在1%以内次数，flag记录连续性，连续5次为true
        // 连续5次差值比例在1%以内，break
        int cnt = 0;
        bool flag = false;
        double pre = 0;

        for (int j=0; j<50; ++j){
            int K = M;
            int N = M;

            CSFloat *A1 = static_cast<CSFloat *>(std::malloc(M*K * sizeof(CSFloat)));

            for(int i=0;i<M*K;++i)
            {
                A1[i] = CSFloat(rand()) / CSFloat(RAND_MAX);
            }
            CSMatrix A(M, K, A1);
            // A.Print();

            CSFloat *A2 = static_cast<CSFloat *>(std::malloc(K*N * sizeof(CSFloat)));

            for(int i=0;i<K*N;++i)
            {
                A2[i] = CSFloat(rand()) / RAND_MAX;
            }
            CSMatrix B(K, N, A2);
            // B.Print();

            CSFloat *A3 = static_cast<CSFloat *>(std::malloc(M*N * sizeof(CSFloat)));

            for(int i=0;i<M*N;++i)
            {
                A3[i] = CSFloat(0);
            }
            CSMatrix C(M, N, A3);
            //C.Print();
            
            int repeat = 20;
            
            Kokkos::Timer timer;
            for(int i=0;i<repeat;++i){
                A.GEMM(B, C);
            }
            Kokkos::fence();
            double time = timer.seconds();
            

            /// \brief 带宽计算公式
            double Gbytes = repeat*1.0e-9*(2.0*M*N*K - M*K)/time;
            
            std::cout<<"CSMatrix A size: "<<setiosflags(ios::left)<<setw(15)<<M*K<<setw(12)<<"CSMatrix B size: "
                        <<setw(15)<<K*N<<"GFlops: "<<Gbytes<<std::endl;

            std::free(A1);
            std::free(A2);
            std::free(A3);

            if(abs(Gbytes-pre)/pre < 0.05){
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
        
            M *= 1.2;
        }
    }
    ChipSum::Common::Finalize();
}