/*
 * @Author: your name
 * @Date: 2021-08-09 13:58:38
 * @LastEditTime: 2021-08-09 15:08:16
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /lky/git/ChipSum/chipsum/numeric/impl/scalar_serial_impl.hpp
 */

#ifndef __CHIPSUM_SCALAR_SERIAL_IMPL_HPP__
#define __CHIPSUM_SCALAR_SERIAL_IMPL_HPP__

#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"
#include "../../backend/backend.hpp"




namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Scalar_Traits<ScalarType,SizeType,ChipSum::Backend::Serial,Props...>:
        public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::Serial>
{
    using scalar_type = ScalarType;
    
    using device_scalar_value_type = ScalarType;

};



namespace Impl{
namespace Scalar{
template<typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void Create(ScalarType&){
    return;
}

template<typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void Create(const ScalarType s,ScalarType& r){
    r=s;
}

template<typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void DeepCopy(const ScalarType s,ScalarType& r){
    r=s;
}



template <typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE ScalarType GetItem(ScalarType s){
    return s;
}


template <typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void GetItem(const ScalarType s,ScalarType& r){
    r=s;
}



}
}

}
}


#endif