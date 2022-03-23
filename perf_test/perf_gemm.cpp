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
        CSFloat tmp;
        tmp = CSFloat(rand()) / CSFloat(RAND_MAX);
        cout<<tmp<<endl;
        cout<<"type: "<<typeid(CSFloat(rand()) / CSFloat(RAND_MAX)).name()<<endl;
        int M = 10;

        for (int j=0; j<50; ++j){
            //M += 10;
            int K = M;
            int N = M;

            CSFloat *A1 = static_cast<CSFloat *>(std::malloc(M*K * sizeof(CSFloat)));
            // Kokkos::View<CSFloat **> tmp = v_unmanged(A1, M, K);

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


            //A.GEMM(B, C);
            //C.Print();

            //(A*B*3).Print();
            
            int repeat = 20;
            /// \brief 暂时用Kokkos的Timer充数吧
            Kokkos::Timer timer;
            for(int i=0;i<repeat;++i){
                A.GEMM(B, C);
            }
            Kokkos::fence();
            double time = timer.seconds();
            

            /// \brief 带宽计算公式
            double Gbytes = repeat*1.0e-9*(2.0*M*N*K - M*K)/time;
            // std::cout<<"---------------------ChipSum Perf Test"
            //     "---------------------"<<std::endl;  
            // std::cout<<M<<std::endl;
            // std::cout<<"Dense matrix GEMM performance : "<<Gbytes<<" GFlops"<<std::endl;
            if(j==0){
                std::cout<<"---------------------ChipSum Perf Test"
                    "---------------------"<<endl
                    <<"CSMatrix size, CSMatrix size, GFlops: "<<std::endl;
            }
            //cout<<i<<endl;
            std::cout<<setiosflags(ios::left)<<setw(15)<<M*K<<setw(12)<<K*N<<Gbytes<<std::endl;
            
            std::free(A1);
            std::free(A2);
            std::free(A3);
            // if (i<30)
            //     M *= 1.2;
            // else
            //     M *= 1.6;
            M *= 1.2;
        }
    }
    ChipSum::Common::Finalize();
}