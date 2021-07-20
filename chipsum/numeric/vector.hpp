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


//    static_assert (
//    std::_or__<>!std::is_same<BackendType,ChipSum::Backend::BuiltinSerial> && , )
    using traits = Vector_Traits<ScalarType,SizeType,BackendType,Props...>;

    using data_type = typename traits::vector_type;
    using const_data_type_reference = typename traits::const_vector_type_reference;

    using size_type = typename traits::size_type;
    using scalar_type = typename traits::scalar_type;
    using scalar_type_reference = typename traits::scalar_type_reference;





private:

    data_type __data;
    size_type __size;

public:

    Vector(){for(int i=0;i<10;++i)__data.push_back(ScalarType(i));}


    inline void setData(const scalar_type* data,const SizeType& size){
//        ChipSum::Numeric::Impl::VectorCopy(data,size,)
    }


    inline const_data_type_reference getData(){return __data;}



    inline void Dot(Vector& v,scalar_type_reference r){



        ChipSum::Numeric::Impl::Dot<ScalarType,SizeType,BackendType>(getData(),v.getData(),r);
    }


};






}
}
#endif // VECTOR_HPP
