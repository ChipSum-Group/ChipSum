#ifndef __VECTOR_SERIAL_IMPL_HPP__
#define __VECTOR_SERIAL_IMPL_HPP__

#include <vector>
#include "../numeric_traits.hpp"


namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUSerialBackend,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUSerialBackend,Props...>
{
    using vector_type = typename std::vector<ScalarType>;
    using vector_type_reference = typename std::vector<ScalarType>::reference ;
    using const_vector_type_reference = typename std::add_const<vector_type_reference>::type;
    using scalar_type = typename std::vector<ScalarType>::value_type;
    using size_type = typename std::vector<ScalarType>::size_type;

};










namespace Impl {






template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>

using traits = Vector_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUSerialBackend,Props...>;

void Dot<ScalarType,SizeType,BackendType,Props...>(typename traits::const_vector_type_reference a,typename traits::const_vector_type_reference b){



}


}

}

}

#endif
