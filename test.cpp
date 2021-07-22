#include <iostream>


#include <vector>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "ChipSumConfig.h"
#include "chipsum/numeric/vector.hpp"
#include "chipsum/numeric/impl/vector_serial_impl.hpp"
#include "chipsum/backend/backend.hpp"


using namespace std;



int main(int argc,char* argv[])
{
    Kokkos::initialize();
    {

    std::cout << " Version " << ChipSum_VERSION_MAJOR << "."
              << ChipSum_VERSION_MINOR << std::endl;


    double* v1 = (double*)malloc(10*sizeof (double));

    double* v2 = (double*)malloc(10*sizeof (double));

    for(int i=0;i<10;++i){
        v1[i]=double(i);
        v2[i]=double(i);
    }

    typedef ChipSum::Numeric::Vector<double,int,ChipSum::Backend::KokkosKernels> Vector;

    Vector a(v1,10);
    Vector b(v2,10);
    double r;
    a.Dot(b,r);

    cout<<r<<endl;


    Kokkos::fence();
    free(v1);
    free(v2);

}
    Kokkos::finalize();

}
