
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

#include "ChipSum.hpp"


#define N 5
#define M 3




int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {
        
        double *v1 = static_cast<double *>(std::malloc(N * sizeof(double)));
        double *v2 = static_cast<double *>(std::malloc(N * sizeof(double)));
        
        for (int i = 0; i < N; ++i) {
            v1[i] = double(i);
            v2[i] = double(i);
        }
        
        Vector a( N,v1); // a = {0,1,2,3,4}
        a.Print();
        Vector b(N,v2); // b = {0,1,2,3,4}
        
        a += b; // a = {0,2,4,6,8}
        
        a += b; // a = {0,3,6,9,12}
        
        a *= 1.0; // a = {0,3,6,9,12}
        
        a.Print();
        
        //    for(std::size_t i=0;i<N;++i) b(i) = 0.0; /* for operator()
        b *= 1.0;  /*  */
        b.Print(); // b = {0.,1.,2., 3. , 4.}
        
        //    a -= a; // if uncomment, all results below turns 0
        
        //    auto c = a+b;
        
        cout << a.Norm1() << endl; // 30
        
        cout << a.Norm2() << endl; // 16.4317



        a.Print();
        a.AXPBY(b,13.0);// a = a+13*b {0 , 16,32,48,64}


        a.Print();

        a.AXPBY(b,1,13); // a=13*a+1*b {0,209,418,627,836}
        a.Print();

        Scalar r ;


        r = 2.7;
        r.Print();
        r = a.Dot(a);
        cout<<r()<<endl; // large, 1.31043e6
        r.Print(); // 1.31043e6

        a.Norm1(r);
        r.Print();



        a.Norm1(r);

        cout<<r()<<endl;

        cout<<a.Norm2()<<endl;

        cout<<a.NormInf()<<endl;


        double *A1 = static_cast<double *>(std::malloc(N*N * sizeof(double)));

        for(int i=0;i<N*N;++i)
        {

            A1[i] = 0;

        }
        for(int i=0;i<N;++i)
        {

            A1[i*N+i] = 1;

        }
        A1[0] = 2;

        double *A2 = static_cast<double *>(std::malloc(N*M * sizeof(double)));

        ;

        Matrix A(N, N,A1);
        A.Print();

        auto y = A*a;
        (A*a).Print();



        for (int i=0;i<N*M;++i) {
            A2[i] = 0;
        }

        for(int i=0;i<M;++i)
        {

            A2[i*M+i] = 2;

        }



        Matrix B(N,M,A2);
        B.Print();

        (A*B*3).Print();


        /*
                *
                *  |  1  0  2  3  0  |
                *  |  0  4  0  5  0  |
                *  |  2  0  6  0  7  |
                *  |  3  5  0  8  0  |
                *  |  0  0  7  0  9  |
                *
                *  nrows = 5
                *  ncols = 5
                *  annz = 13
                *  row_map = {0,3,5,8,11,13}
                *  col_map = {0,2,3,1,3,0,2,4,0,1,3,2,4}
                *  values = {1,2,3,4,5,2,6,7,3,5,8,7,9}
                *
                *  x = {1,1,1,1,1}
                *  y = {6,9,15,16,16}
                *
               */
        CSInt m = 5;
        CSInt n = 5;
        
        CSInt nrows = m;
        CSInt ncols = n;
        CSInt annz = 13;
        CSInt *row_map = (CSInt *)malloc(6 * sizeof(CSInt));
        CSInt *col_map = (CSInt *)malloc(13 * sizeof(CSInt));
        double *values = (double *)malloc(13 * sizeof(double));
        row_map[0] = 0;
        row_map[1] = 3;
        row_map[2] = 5;
        row_map[3] = 8;
        row_map[4] = 11;
        row_map[5] = 13;
        col_map[0] = col_map[5] = col_map[8] = 0;
        col_map[1] = col_map[6] = col_map[11] = 2;
        col_map[2] = col_map[4] = col_map[10] = 3;
        col_map[3] = col_map[9] = 1;
        col_map[7] = col_map[12] = 4;
        
        for (int i = 0; i < 5; ++i)
            values[i] = double(i + 1);
        
        values[5] = 2;
        values[6] = 6;
        values[7] = 7;
        values[8] = 3;
        values[9] = 5;
        values[10] = 8;
        values[11] = 7;
        values[12] = 9;

        CSR s_A(nrows,ncols,annz,row_map,col_map,values);

        b.Print();
        (s_A*b).Print();


        std::free(row_map);
        std::free(col_map);
        std::free(values);
        
        //        const CSInt m = 500;
        //        const CSInt n = 500;
        //        CSInt nrows = m;
        //        CSInt ncols = n;
        //        CSInt annz = m;
        //        CSInt *row_map = (CSInt *)malloc((m + 1) * sizeof(CSInt));
        //        CSInt *col_map = (CSInt *)malloc(std::min(m,n) * sizeof(CSInt));
        //        double *values = (double *)malloc(std::min(m,n) * sizeof(double));
        
        //        for (CSInt i = 0; i < std::min(m,n); ++i) {

        //            col_map[i] = i;
        //            values[i] = 1;
        //        }
        //        for (CSInt i = 0; i <= m; ++i) {

        //            row_map[i] = i;
        //        }
        

        
        //        bb.Print(); // {13,15,40,24,50}
        
        //        Kokkos::fence();
        
        //        cout<<endl;
        
        //        free(row_map);
        //        free(col_map);
        //        free(values);
        
        //        //    for(size_t i=0;i<5;++i){ /* for operator() test */
        //        //        mat(i,i) = 1.0;
        //        //    }
        
        //        double* mat_data = static_cast<double*>(std::malloc(25*sizeof
        //        (double)));
        
        //        for(int i=0;i<25;++i){
        //            mat_data[i] = 0.0;
        //        }
        
        //        for(int i=0;i<5;++i){
        //            mat_data[5*i+i] = 1.0;
        //        }
        
        //        mat_data[0] = 2.0;
        
        //        Matrix mat(5,5,mat_data);
        
        //        mat.Print();
        
        //        /*     ↑
        //        *
        //        *  |  2  0  0  0  0  |
        //        *  |  0  1  0  0  0  |
        //        *  |  0  0  1  0  0  |
        //        *  |  0  0  0  1  0  |
        //        *  |  0  0  0  0  1  |
        //        */
        
        //        (mat*bb).Print(); // {13,15,40,24,50}
        
        //        (B*mat).Print();
        
        //        /*     ↑
        //        *
        //        *  |  2  0  2  3  0  |
        //        *  |  0  4  0  5  0  |
        //        *  |  4  0  6  0  7  |
        //        *  |  6  5  0  8  0  |
        //        *  |  0  0  7  0  9  |
        //        */
        
        //        std::free(mat_data);
        
        std::free(v1);
        std::free(v2);
        std::free(A1);
        std::free(A2);
    }
    ChipSum::Common::Finalize();
}
