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

    ChipSum::Numeric::Vector<double,int,ChipSum::Backend::KokkosKernels> a;
    ChipSum::Numeric::Vector<double,int,ChipSum::Backend::KokkosKernels> b;
    double r;
    a.Dot(b,r);

    cout<<r<<endl;
}
    Kokkos::finalize();

}
