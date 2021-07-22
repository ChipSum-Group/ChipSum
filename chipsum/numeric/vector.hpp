#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <vector>
#include <type_traits>

//#include <iostream>
//using namespace std;

#include "numeric_traits.hpp"
#include "impl/vector_serial_impl.hpp"
#include "impl/vector_kokkoskernels_impl.hpp"

#include <KokkosBlas1_fill.hpp>

static int vector_name_counter;

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



    explicit inline Vector(scalar_type* data,const SizeType& size){


        ChipSum::Numeric::Impl::Vector::CreateVector<ScalarType,SizeType>(size,__data);
        ChipSum::Numeric::Impl::Vector::FillVector<ScalarType,SizeType>(data,size,__data);

    }


    inline void setData(const scalar_type* data,const SizeType& size){
        ChipSum::Numeric::Impl::Vector::FillVector(data,size,__data);
    }


    inline const_data_type_reference getData(){return __data;}



    inline void Dot(Vector& v,scalar_type_reference r){
        ChipSum::Numeric::Impl::Vector::Dot<
                ScalarType,SizeType,BackendType,Props...
                >(getData(),v.getData(),__size,r);
    }


};

//int ChipSum::Numeric::Vector::__name_counter = 0;




}
}
#endif // VECTOR_HPP
