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


template<typename ...Props>
struct S1{

    using t0 = int;

};

template<typename T,typename ...Props>
struct S1<T,Props...>{

    using t1 = double;
};




int main(int argc,char* argv[])
{

    Kokkos::initialize();
    {

    double* v1 = (double*)malloc(10*sizeof (double));

    double* v2 = (double*)malloc(10*sizeof (double));

    for(int i=0;i<10;++i){
        v1[i]=double(i);
        v2[i]=double(i);
    }


    S1<int,double>::t1 a1 = 999;

    cout<<a1<<endl;


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
}
    Kokkos::finalize();

}
