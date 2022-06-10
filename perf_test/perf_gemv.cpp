/* * * * * * * * * * * * * * * * * * * * *
 *   File:     test.cpp
 *   Author:   Li Kunyun
 *   group:    CDCS-HPC
 *   Time:     2021-07-28
 * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using namespace std;

#include <type_traits>
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

        for (int j=0; j<50; j++){
            //M += 10;
            int N = M;

            CSFloat *A1 = static_cast<CSFloat *>(std::malloc(M*N * sizeof(CSFloat)));

            for(int i=0;i<M*N;++i)
            {
                A1[i] = CSFloat(rand()) / RAND_MAX;
            }
            CSMatrix A(M, N, A1);

            CSFloat *A2 = static_cast<CSFloat *>(std::malloc(N * sizeof(CSFloat)));

            for(int i=0;i<N;++i)
            {
                A2[i] = CSFloat(rand()) / RAND_MAX;
            }
            CSVector x(N, A2);

            CSFloat *A3 = static_cast<CSFloat *>(std::malloc(M * sizeof(CSFloat)));

            for(int i=0;i<M;++i)
            {
                A3[i] = CSFloat(0);
            }
            CSVector r(M, A3);

            
            int repeat = 20;
            /// \brief 暂时用Kokkos的Timer充数吧
            Kokkos::Timer timer;
            for(int i=0;i<repeat;++i){
                A.GEMV(x, r);
            }
            Kokkos::fence();
            double time = timer.seconds();
            

            /// \brief 带宽计算公式
            double Gbytes = repeat*1.0e-9*(2.0*M*N - M)/time;
            
            cout<<"CSMatrix size: "<<setiosflags(ios::left)<<setw(15)<<M*M<<setw(12)<<"CSVector size: "
                        <<setw(12)<<N<<"GFlops: "<<Gbytes<<endl;
            
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