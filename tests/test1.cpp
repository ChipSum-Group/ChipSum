
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

#include <Kokkos_Core.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include "../chipsum/common/coo_reader.hpp"

#define N 5
#define M 3


struct ParallelPrintf_matrix {
    
    CSMatrix _data;

    ParallelPrintf_matrix(CSMatrix data) : _data(data){}

    CHIPSUM_SPECIAL_INLINE
    void operator() (const int64_t i, const int64_t j) const {
        printf("%f** ", _data.Item(i,j)); 
        // printf("%f** ", _data[i][j]); 
    }
};

struct ParallelPrintf_vec {

    ParallelPrintf_vec(double *v):_vec(N,v){   
        _vec.Print();
    }

    CHIPSUM_SPECIAL_INLINE
    void operator ()(const int i ) const {
        _vec[i] = 12;
        printf("AAA: %f\n", _vec.Item(i));

    }

    private:
    CSVector _vec;
};


int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {
        
        double *v1 = static_cast<double *>(std::malloc(N * sizeof(double)));
        double *v2 = static_cast<double *>(std::malloc(N * sizeof(double)));
        
        for (int i = 0; i < N; ++i) {
            v1[i] = double(i);
            v2[i] = double(i);
        }
        
        CSVector a( N,v1); // a = {0,1,2,3,4}
        a.Print();
        CSVector b(N,v2); // b = {0,1,2,3,4}
        
        // a.HostToDevice();
        // a.DeviceToHost();
        
        a += b; // a = {0,2,4,6,8}

        cout<<"vector host index "<<a(1)<<std::endl;
        Kokkos::parallel_for("dev_p",N, ParallelPrintf_vec(v1));
        
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

            A1[i] = -1;

        }
        for(int i=0;i<N;++i)
        {

            A1[i*N+i] = 1;

        }
        A1[0] = 2;

        double *A2 = static_cast<double *>(std::malloc(N*M * sizeof(double)));

        CSMatrix A(N, N, A1);
        A.Print();

        A.Norm();
        A.Print();
        // A.LeakyRelu();
        A.Softmax();

        A.Print();

        using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        ParallelPrintf_matrix functor(A);
        Kokkos::parallel_for( "device_out", mdrange_policy({0,0}, {5,5}), functor);

        A(0,0)=0;
        A.HostToDevice();
        A.DeviceToHost();
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



        CSMatrix B(N,M,A2);
        B.Print();

        double *A3 = static_cast<double *>(std::malloc(N*M * sizeof(double)));
        for (int i=0;i<N*M;++i) {
            A3[i] = 0;
        }
        CSMatrix C(N,M,A3);

        double *A4 = static_cast<double *>(std::malloc(M * sizeof(double)));
        for (int i=0;i<M;++i) {
            A4[i] = -100;
        }
        CSVector bias(M, A4);

        A.Dense("N", "N", B, C);
        A.Dense("N", "N", B, bias, C);
        C.Print();
        cout<<"dense result"<<endl;
        
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
        CSFloat *values = (double *)malloc(13 * sizeof(double));
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
/*
        CSInt nv = 0, ne = 0;
          CSInt *xadj, *adj;
          CSFloat *ew;

        std::string mtx_file("/home/lky/bin/1.mtx");
        KokkosKernels::Impl::read_matrix<CSInt,CSInt,CSFloat>(&nv,&ne,&xadj,&adj,&ew,mtx_file.data());

        CSR s_M(nv,nv,ne,xadj,adj,ew);

        s_M.SavePatternFig("mtx.PNG");
*/
        CSInt *row_map2 = (CSInt *)malloc(5 * sizeof(CSInt));
        CSInt *col_map2 = (CSInt *)malloc(5 * sizeof(CSInt));
        CSFloat *values2 = (double *)malloc(5 * sizeof(double));
        row_map2[0] = 0;
        row_map2[1] = 1;
        row_map2[2] = 2;
        row_map2[3] = 3;
        row_map2[4] = 4;

        col_map2[0] = 0;
        col_map2[1] = 1;
        col_map2[2] = 2;
        col_map2[3] = 3;
        col_map2[4] = 4;

        values2[0] = 1;
        values2[1] = 2;
        values2[2] = 3;
        values2[3] = 4;
        values2[4] = 5;
        COO s_C(nrows,ncols,5,row_map2,col_map2,values2);
        s_C.Insert(1,2,7.0);
        std::vector<int> t_row_map, t_col_map;
        std::vector<CSFloat> t_values;
        // s_C.GetCrsData(t_row_map, t_col_map, t_values);


        //read coo
        /* int nr, nc, nz;
        int *rm, *cm;
        CSFloat *vals;
        std::vector<int> t_row_map1, t_col_map1;
        std::vector<CSFloat> t_values1;
        ChipSum::Common::coo_reader(nr, nc, nz, rm, cm, vals, "../../data/A.mtx" );
        COO sparse_coo(nr,nc,nz,rm,cm,vals);
        //coo to csr
        sparse_coo.GetCrsData(t_row_map1, t_col_map1, t_values1);
        CSR sparse_csr(nr,nc,nz,t_row_map1.data(),t_col_map1.data(),t_values1.data()); */

        


        // delete rm;
        // delete cm;
        // delete vals;

        std::free(row_map);
        std::free(col_map);
        std::free(values);
        std::free(row_map2);
        std::free(col_map2);
        std::free(values2);
        

        std::free(v1);
        std::free(v2);
        std::free(A1);
        std::free(A2);
        std::free(A3);
        std::free(A4);
    }
    ChipSum::Common::Finalize();
}