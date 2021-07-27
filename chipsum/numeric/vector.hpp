#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <vector>
#include <type_traits>

//#include <iostream>
//using namespace std;

#include "numeric_traits.hpp"
#include "impl/vector_serial_impl.hpp"
#include "impl/vector_kokkoskernels_impl.hpp"


static int vector_name_counter;

namespace ChipSum {
namespace Numeric {

template<typename ...Props>
class Vector;

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
class Vector<ScalarType,SizeType,BackendType,Props...>{

public:


    using traits = Vector_Traits<ScalarType,SizeType,Props...>;
    using vtraits = Vector_Traits<ScalarType,SizeType,BackendType,Props...>;

    using data_type = typename vtraits::vector_type;
    using data_type_reference = typename std::add_lvalue_reference<data_type>::type ;
    using const_data_type_reference = typename std::add_const<data_type_reference>::type;

    // Not core
    using size_type = typename std::remove_const<typename vtraits::size_type>::type;
    using const_size_type = typename std::add_const<size_type>::type;
    using size_type_reference = typename std::add_lvalue_reference<size_type>::type;
    using const_size_type_reference = typename add_const<size_type_reference>::type;

    using scalar_type = typename traits::nonconst_scalar_type;
    using scalar_type_reference = typename traits::nonconst_scalar_type_reference;



private:

    data_type __data;
    size_type __size;

private:

    //    void privateSample(){}

protected:

    //    void _protectedSample(){}

public:

    //    void PublicSample(){}

    explicit CHIPSUM_FUNCTION_INLINE Vector(const data_type& data,const size_type size):
        __data(data),__size(size){}


    explicit CHIPSUM_FUNCTION_INLINE Vector(scalar_type* data,const SizeType& size)
        :__size(size)
    {
        ChipSum::Numeric::Impl::Vector::Create<ScalarType,SizeType>(size,__data);
        ChipSum::Numeric::Impl::Vector::Fill(data,size,__data);
    }


    CHIPSUM_FUNCTION_INLINE void SetData(const scalar_type* data,const SizeType& size){
        ChipSum::Numeric::Impl::Vector::Fill(data,size,__data);
    }


    //    CHIPSUM_FUNCTION_INLINE const_size_type

    CHIPSUM_FUNCTION_INLINE const_data_type_reference GetData(){return __data;}

    CHIPSUM_FUNCTION_INLINE const_size_type_reference GetSize(){return __size;}



    CHIPSUM_FUNCTION_INLINE void Dot(Vector& v,scalar_type_reference r){
        ChipSum::Numeric::Impl::Vector::
                Dot(GetData(),v.GetData(),__size,r);
    }


    CHIPSUM_FUNCTION_INLINE Vector operator*(ScalarType s){

        data_type out_data;

        ChipSum::Numeric::Impl::Vector::Create(__size,out_data);

        ChipSum::Numeric::Impl::Vector::Scal<ScalarType,SizeType,Props...>(out_data,s,GetData());

        Vector out(out_data,__size);

        return out;
    }

    CHIPSUM_FUNCTION_INLINE Vector& operator*=(ScalarType s){
        ChipSum::Numeric::Impl::Vector::Scal<ScalarType,SizeType,Props...>(__data,s,GetData());
        return *this;
    }


    CHIPSUM_FUNCTION_INLINE ScalarType Norm1(){
        return ChipSum::Numeric::Impl::Vector::Norm1<ScalarType,SizeType,Props...>(__data);
    }


    CHIPSUM_FUNCTION_INLINE ScalarType Norm2(){
        return ChipSum::Numeric::Impl::Vector::Norm2<ScalarType,SizeType,Props...>(__data);
    }




};



template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
CHIPSUM_FUNCTION_INLINE Vector<ScalarType,SizeType,BackendType,Props...> operator*(ScalarType s,Vector<ScalarType,SizeType,BackendType,Props...>& v){
    return v*s;
}



} // End namespace Numeric
} // End namespace ChipSum



#endif // VECTOR_HPP
