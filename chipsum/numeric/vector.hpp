#ifndef __CHIPSUM_NUMERIC_VECTOR_HPP__
#define __CHIPSUM_NUMERIC_VECTOR_HPP__

#include <vector>
#include <type_traits>

/* To debug */
//#include <iostream>
//using namespace std;

#include "numeric_traits.hpp"
#include "impl/vector_serial_impl.hpp"
#include "impl/vector_kokkoskernels_impl.hpp"




namespace ChipSum {
namespace Numeric {

template<typename ...Props>
class Vector;

template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
class Vector<ScalarType,SizeType,BackendType,Props...>{

public:


    using traits = Vector_Traits<ScalarType,SizeType,BackendType,Props...>;


    using vector_type = typename traits::vector_type;
    using size_type = typename traits::size_type;
    using size_type_reference = typename std::add_lvalue_reference<size_type>::type;
    using const_size_type_reference = typename std::add_const<size_type_reference>::type;

    using vector_type_reference = typename std::add_lvalue_reference<vector_type>::type;
    using const_vector_type_reference = typename std::add_const<vector_type_reference>::type;

private:

    vector_type __data;
    size_type __size;

private:

    //    void privateSample(){}

protected:

    //    void _protectedSample(){}

public:

    //    void PublicSample(){}

    /**
     * @brief Vector
     * @param data
     * @param size
     */
    explicit CHIPSUM_FUNCTION_INLINE Vector(const vector_type& data,const size_type size):
        __data(data),__size(size){}


    /**
     * @brief Vector
     * @param data
     * @param size
     */
    explicit CHIPSUM_FUNCTION_INLINE Vector(typename traits::nonconst_scalar_type* data,const SizeType& size)
        :__size(size)
    {
        ChipSum::Numeric::Impl::Vector::Create<ScalarType,SizeType>(size,__data);
        ChipSum::Numeric::Impl::Vector::Fill(data,size,__data);
    }


    /**
     * @brief SetData
     * @param data
     * @param size
     */
    CHIPSUM_FUNCTION_INLINE void SetData(typename traits::nonconst_scalar_type* data,const SizeType& size){
        ChipSum::Numeric::Impl::Vector::Fill(data,size,__data);
    }


    /**
     * @brief GetData
     * @return
     */

    CHIPSUM_FUNCTION_INLINE const_vector_type_reference GetData(){return __data;}

    /**
     * @brief GetSize
     * @return
     */
    CHIPSUM_FUNCTION_INLINE const_size_type_reference GetSize(){return __size;}


    /**
     * @brief Dot
     * @param v
     * @param r
     */

    CHIPSUM_FUNCTION_INLINE void Dot(Vector& v,typename traits::nonconst_scalar_type_reference r){
        ChipSum::Numeric::Impl::Vector::
                Dot(GetData(),v.GetData(),__size,r);
    }


    /**
     * @brief operator *
     * @param s
     * @return
     */
    CHIPSUM_FUNCTION_INLINE Vector operator*(ScalarType s){

        vector_type out_data;

        ChipSum::Numeric::Impl::Vector::Create(__size,out_data);

        ChipSum::Numeric::Impl::Vector::Scal<ScalarType,SizeType,Props...>(out_data,s,GetData());

        Vector out(out_data,__size);

        return out;
    }

    /**
     * @brief operator *=
     * @param s
     * @return
     */
    CHIPSUM_FUNCTION_INLINE Vector& operator*=(ScalarType s){
        ChipSum::Numeric::Impl::Vector::Scal<ScalarType,SizeType,Props...>(__data,s,GetData());
        return *this;
    }

    /**
     * @brief operator +
     * @param s
     * @return
     */
    CHIPSUM_FUNCTION_INLINE Vector& operator+(ScalarType s){
        //TODO
    }

    /**
     * @brief operator +=
     * @param s
     * @return
     */
    CHIPSUM_FUNCTION_INLINE Vector& operator+=(ScalarType s){
        //TODO
    }

    /**
     * @brief operator -
     * @param s
     * @return
     */
    CHIPSUM_FUNCTION_INLINE Vector& operator-(ScalarType s){
        //TODO
    }

    /**
     * @brief operator -=
     * @param s
     * @return
     */
    CHIPSUM_FUNCTION_INLINE Vector& operator-=(ScalarType s){
        //TODO
    }



    /**
     * @brief Norm1
     * @return
     */
    CHIPSUM_FUNCTION_INLINE ScalarType Norm1(){
        return ChipSum::Numeric::Impl::Vector::Norm1<ScalarType,SizeType,Props...>(__data);
    }


    /**
     * @brief Norm2
     * @return
     */
    CHIPSUM_FUNCTION_INLINE ScalarType Norm2(){
        return ChipSum::Numeric::Impl::Vector::Norm2<ScalarType,SizeType,Props...>(__data);
    }


};



/**
 *
 */
template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
CHIPSUM_FUNCTION_INLINE Vector<ScalarType,SizeType,BackendType,Props...>
operator*(ScalarType s,Vector<ScalarType,SizeType,BackendType,Props...>& v){
    return v*s;
}



} // End namespace Numeric
} // End namespace ChipSum



#endif // VECTOR_HPP
