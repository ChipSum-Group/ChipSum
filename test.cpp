#include <iostream>
#include <Testconfig.h>

#include <vector>
#include <type_traits>

#include "chipsum/numeric/vector.hpp"
#include "chipsum/numeric/impl/vector_serial_impl.hpp"

using namespace std;

template<typename ...>
void dot(...);

template<typename T>
void dot(T& i,T& j){

    cout<<sizeof(T)<<endl;
}

int main(int argc,char* argv[])
{

    std::cout << " Version " << ChipSum_VERSION_MAJOR << "."
              << ChipSum_VERSION_MINOR << std::endl;


    //    ChipSum::Numeric::Operator_Traits<double,int>::const_size_type a=1.2;

    ChipSum::Numeric::Vector_Traits<double,int,ChipSum::Backend::KokkosKernels> b;


    ChipSum::Numeric::Vector_Traits<double,int,ChipSum::Backend::DefaultBackend> v1;
    ChipSum::Numeric::Vector_Traits<double,int,ChipSum::Backend::DefaultBackend> v2;


//    ChipSum::Numeric::Impl::Dot();
    double a;
//    ChipSum::Numeric::Dot<double,int>(v1,v2,a);


    //    cout<<"IS BASE OF:"<<std::is_base_of<A<>,B<>>::value<<endl;


    std::vector<int> v;
}
