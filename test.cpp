#include <iostream>


#include <vector>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "ChipSumConfig.h"
#include "chipsum/numeric/vector.hpp"
#include "chipsum/numeric/impl/vector_serial_impl.hpp"
#include "chipsum/backend/backend.hpp"


using namespace std;
typedef ChipSum::Numeric::Vector<double,int,ChipSum::Backend::KokkosKernels> Vector;


template<typename T1,typename T2,typename ...Props>
struct S1{

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
