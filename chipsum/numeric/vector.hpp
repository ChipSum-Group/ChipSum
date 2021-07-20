#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__

#include <vector>
#include <type_traits>

//#include <iostream>
//using namespace std;
#include "numeric_traits.hpp"
#include "impl/vector_serial_impl.hpp"



namespace ChipSum {
namespace Numeric {






template<typename ...Props>
class Vector;

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
class Vector<ScalarType,SizeType,BackendType,Props...>{

    using traits = Vector_Traits<ScalarType,SizeType,BackendType,Props...>;

    using data_type = typename traits::vector_type;
    using size_type = typename traits::size_type;



private:

     data_type __data;
     size_type __size;

public:

     void Dot(Vector& v){

     }


};





}
}
#endif // VECTOR_HPP
