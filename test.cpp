#include <iostream>
#include <Testconfig.h>

#include <vector>
#include <type_traits>

#include "chipsum/numeric/vector.hpp"
#include "chipsum/numeric/impl/vector_serial_impl.hpp"
#include "chipsum/backend/backend.hpp"

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

    ChipSum::Numeric::Vector<double,int,ChipSum::Backend::BuiltinSerial> a;
    ChipSum::Numeric::Vector<double,int,ChipSum::Backend::BuiltinSerial> b;
    double r;
    a.Dot(b,r);

    cout<<r<<endl;


}
