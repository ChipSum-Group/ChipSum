#include <iostream>


#include <vector>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "ChipSumConfig.h"
#include "chipsum/numeric/vector.hpp"
#include "chipsum/numeric/impl/vector_serial_impl.hpp"
#include "chipsum/backend/backend.hpp"
#include "chipsum/numeric/sparse_matrix.hpp"


using namespace std;
typedef ChipSum::Numeric::Vector<double,int,ChipSum::Backend::KokkosKernels> Vector;





int main(int argc,char* argv[])
{

    Kokkos::initialize();
    {

    double* v1 = (double*)malloc(10*sizeof (double));

    double* v2 = (double*)malloc(10*sizeof (double));

    for(int i=0;i<10;++i){
        v1[i]=1;
        v2[i]=1;
    }





    Vector a(v1,10);
    Vector b(v2,10);


    a = 1.0*a;


    cout<<a.Norm1()<<endl;

    cout<<a.Norm2()<<endl;
    double r;
    a.Dot(a,r);

    Kokkos::fence();
    cout<<r<<endl;

    free(v1);
    free(v2);



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
     */

    size_t nrows = 5;
    size_t ncols = 5;
    size_t annz = 13;
    size_t* row_map = (size_t*)malloc(6*sizeof (size_t));
    size_t* col_map = (size_t*)malloc(13*sizeof (size_t));
    double* values = (double*)malloc(13*sizeof (double));
    row_map[0]=0;row_map[1]=3;row_map[2]=5;row_map[3]=8;row_map[4]=11;row_map[5]=13;
    col_map[0]=col_map[5]=col_map[8] = 0;
    col_map[1]=col_map[6]=col_map[11] = 2;
    col_map[2]=col_map[4]=col_map[10] = 3;
    col_map[3]=col_map[9] = 1;
    col_map[7]=col_map[12] = 4;

    for(int i=0;i<5;++i) values[i]=double(i+1);

    values[5]=2;values[6]=6;values[7]=7;
    values[8]=3;values[9]=5;values[10]=8;
    values[11]=7;values[12]=9;


    KokkosSparse::CrsMatrix<double,size_t ,default_device> A;
    ChipSum::Numeric::Impl::Sparse::Fill(A,nrows,ncols,annz,row_map,col_map,values);

//    CHIPSUM_FUNCTION_INLINE void Fill(
//            KokkosSparse::CrsMatrix<ScalarType,SizeType,default_device>& A,
//            const SizeType nrows,
//            const SizeType ncols,
//            const SizeType annz,
//            SizeType* row_map,
//            SizeType* col_map,
//            ScalarType* values

//            )

    ChipSum::Numeric::SparseMatrix<double,size_t,ChipSum::Numeric::Csr,ChipSum::Backend::KokkosKernels>
            B(nrows,ncols,annz,row_map,col_map,values);

    Kokkos::View<double*> x("x",5);

    KokkosBlas::fill(x,1.0);

    Kokkos::View<double*> y("y",5);
    ChipSum::Numeric::Impl::Sparse::Spmv(A,x,y);

    // should be {6,9,15,16,16}
    Kokkos::parallel_for(5,KOKKOS_LAMBDA(int i){
                            printf("%.0f,",y(i));
                         });



    printf("%f\n",y(3));

    Kokkos::fence();

    cout<<endl;

    free(row_map);
    free(col_map);
    free(values);


}
    Kokkos::finalize();

}
