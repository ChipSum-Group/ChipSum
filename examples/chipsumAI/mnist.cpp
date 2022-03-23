
/* * * * * * * * * * * * * * * * * * * * *
 *   File:     mnist.cpp
 *   Author:   Xiang Yukai
 *   group:    CDCS-HPC
 *   Time:     2022-02-15
 * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using namespace std;

#include <type_traits>
#include <vector>

#include "../ChipSum.hpp"

#include <Kokkos_Core.hpp>
#include <KokkosKernels_IOUtils.hpp>

#define M 784
#define N 512
#define K 10

void init_input(int index, double *pic, char *pic_28){
    
    // init input
    ifstream pic_in("../../../data/mnist/data/pic_" + to_string(index) + ".txt");
    for(int i = 0; i < M; ++i){
        pic_in >> pic[i];
        pic[i] /= 255;
        if(pic[i]>0){
            pic_28[i] = '#';
        }
        else{
            pic_28[i] = ' ';
        }
    }
    pic_in.close();
    
    // 可视化输入图片
    cout << "******input is****** : "<< index << endl;
    for (std::size_t i = 0; i < 28; ++i) {
        cout << " "
            << "[";
        for (std::size_t j = 0; j < 28 - 1; ++j) {
            cout << pic_28[i*28+j] << " ";
        }
        cout << pic_28[i*28+27] << "]" << endl;
    }
}

void init_ModelParameter(double *w1, double *w2, double *w3, double *b1, double *b2, double *b3){
    
    // init dense1
    ifstream w1_in("../../../data/mnist/parameters/weight_Dense_1.txt");
    for(int i = 0; i < M; ++i){
        for(int j = 0; j < N; ++j){
            w1_in >> w1[i*N+j];
        }
    }
    w1_in.close();

    ifstream b1_in("../../../data/mnist/parameters/bias_Dense_1.txt");
    for(int i = 0; i < N; ++i){
        b1_in >> b1[i];
    }
    b1_in.close();

    // init dense2
    ifstream w2_in("../../../data/mnist/parameters/weight_Dense_2.txt");
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            w2_in >> w2[i*N+j];
        }
    }
    w2_in.close();

    ifstream b2_in("../../../data/mnist/parameters/bias_Dense_2.txt");
    for(int i = 0; i < N; ++i){
        b2_in >> b2[i];
    }
    b2_in.close();

    // init dense3
    ifstream w3_in("../../../data/mnist/parameters/weight_Dense_3.txt");
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < K; ++j){
            w3_in >> w3[i*K+j];
        }
    }
    w3_in.close();

    ifstream b3_in("../../../data/mnist/parameters/bias_Dense_3.txt");
    for(int i = 0; i < K; ++i){
        b3_in >> b3[i];
    }
    b3_in.close();
}


int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {
        double *w1 = static_cast<double *>(std::malloc(M*N * sizeof(double)));
        double *b1 = static_cast<double *>(std::malloc(N * sizeof(double)));

        double *w2 = static_cast<double *>(std::malloc(N*N * sizeof(double)));
        double *b2 = static_cast<double *>(std::malloc(N * sizeof(double)));

        double *w3 = static_cast<double *>(std::malloc(N*K * sizeof(double)));
        double *b3 = static_cast<double *>(std::malloc(K * sizeof(double)));

        double *pic = static_cast<double *>(std::malloc(M * sizeof(double)));
        char *pic_28 = static_cast<char *>(std::malloc(M * sizeof(char)));

        // init model parameters
        init_ModelParameter(w1, w2, w3, b1, b2, b3);

        // init parameters CSMatrix
        CSMatrix weight_Dense_1(M, N, w1);
        CSVector bias_Dense_1(N, b1);

        CSMatrix weight_Dense_2(N, N, w2);
        CSVector bias_Dense_2(N, b2);

        CSMatrix weight_Dense_3(N, K, w3);
        CSVector bias_Dense_3(K, b3);

        // temp variable
        CSMatrix C1(1,N);
        CSMatrix C2(1,N);
        CSMatrix C3(1,K);


        // compute
        for (int index=0; index<10; ++index)
        {
            // get input
            init_input(index, pic, pic_28);
            CSMatrix input(1, M, pic);

            // model
            input.Dense("N", "N", weight_Dense_1, bias_Dense_1, C1);
            C1.Relu();

            C1.Dense("N", "N", weight_Dense_2, bias_Dense_2, C2);
            C2.Relu();

            C2.Dense("N", "N", weight_Dense_3, bias_Dense_3, C3);
            C3.Softmax();

            C3.Argmax();
        }

        
        std::free(w1);
        std::free(w2);
        std::free(w3);
        std::free(b1);
        std::free(b2);
        std::free(b3);
        std::free(pic);
        std::free(pic_28);

    }
    ChipSum::Common::Finalize();
}
