/* * * * * * * * * * * * * * * * * * * * *
*   File:     test.cpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using namespace std;

#include <vector>
#include <type_traits>

#include "ChipSumConfig.h"
#include "chipsum/numeric/vector.hpp"
#include "chipsum/numeric/impl/vector_serial_impl.hpp"
#include "chipsum/backend/backend.hpp"
#include "chipsum/numeric/sparse_matrix.hpp"
#include "chipsum/numeric/dense_matrix.hpp"
#include "chipsum/numeric/scalar.hpp"




typedef ChipSum::Numeric::SparseMatrix<double,size_t,ChipSum::Numeric::SparseTypes::Csr,ChipSum::Backend::DefaultBackend>  Csrm;

int main(int argc,char* argv[])
{

    Kokkos::initialize();
    {




        double* v1 = (double*)malloc(10*sizeof (double));

        double* v2 = (double*)malloc(10*sizeof (double));

        for(int i=0;i<10;++i){
            v1[i]=i;
            v2[i]=i;
        }

        Vector x(v1,10);
        x.Print();

        Vector y(v2,10);

        Scalar s = 3.3;


        free(v1);free(v2);

        cout<<s()<<endl;

        (s*x).Print();

        s = x.Dot(y);

        s.Print();

        Scalar s1;
        x.Dot(y,s1);
        s1.Print();

        cout<<sizeof(Scalar)<<endl;

//        Vector a(v1,10); //a = {0,1,2,3,4,5,6,7,8,9}
//        Vector b(v2,10); //b = {0,1,2,3,4,5,6,7,8,9}


//        a+=b; //a = {0,2,4,6,8,10,12,14,16,18}

//        a = a+b; //a = {0,3,6,9,12,15,18,21,24,27}


//        a = 1.0*a; //a = {0,3,6,9,12,15,18,21,24,27}

//        a.Print();



//        //    for(std::size_t i=0;i<10;++i) b(i) = 0.0; /* for operator() test */
//        b *= 0.0; /*  */
//        b.Print(); // b = {0.,0.,0., ... ,0.}


//        //    a -= a; // if uncomment, all results above turns 0

//        //    auto c = a+b;

//        cout<<a.Norm1()<<endl; // 135

//        cout<<a.Norm2()<<endl; // 50.6458
//        double r;
//        a.Dot(a,r);

//        cout<<r<<endl; // r = 2565

//        free(v1);
//        free(v2);




//        /*
//        *
//        *  |  1  0  2  3  0  |
//        *  |  0  4  0  5  0  |
//        *  |  2  0  6  0  7  |
//        *  |  3  5  0  8  0  |
//        *  |  0  0  7  0  9  |
//        *
//        *  nrows = 5
//        *  ncols = 5
//        *  annz = 13
//        *  row_map = {0,3,5,8,11,13}
//        *  col_map = {0,2,3,1,3,0,2,4,0,1,3,2,4}
//        *  values = {1,2,3,4,5,2,6,7,3,5,8,7,9}
//        *
//        *  x = {1,1,1,1,1}
//        *  y = {6,9,15,16,16}
//        *
//       */

//        size_t nrows = 5;
//        size_t ncols = 5;
//        size_t annz = 13;
//        size_t* row_map = (size_t*)malloc(6*sizeof (size_t));
//        size_t* col_map = (size_t*)malloc(13*sizeof (size_t));
//        double* values = (double*)malloc(13*sizeof (double));
//        row_map[0]=0;row_map[1]=3;row_map[2]=5;row_map[3]=8;row_map[4]=11;row_map[5]=13;
//        col_map[0]=col_map[5]=col_map[8] = 0;
//        col_map[1]=col_map[6]=col_map[11] = 2;
//        col_map[2]=col_map[4]=col_map[10] = 3;
//        col_map[3]=col_map[9] = 1;
//        col_map[7]=col_map[12] = 4;

//        for(int i=0;i<5;++i) values[i]=double(i+1);

//        values[5]=2;values[6]=6;values[7]=7;
//        values[8]=3;values[9]=5;values[10]=8;
//        values[11]=7;values[12]=9;

//        Csrm B(nrows,ncols,annz,row_map,col_map,values);
////        B.GSSmooth();


//        Vector xb(v1,5);

//        Vector bb = B*xb;

//        bb.Print(); // {13,15,40,24,50}

//        Kokkos::fence();

//        cout<<endl;

//        free(row_map);
//        free(col_map);
//        free(values);




//        //    for(size_t i=0;i<5;++i){ /* for operator() test */
//        //        mat(i,i) = 1.0;
//        //    }

//        double* mat_data = static_cast<double*>(std::malloc(25*sizeof (double)));

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

    }
    Kokkos::finalize();

}
