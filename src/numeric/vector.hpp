#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <type_traits>

#include <iostream>
using namespace std;
#include "operator.hpp"



namespace ChipSum {
namespace Numeric {


template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
struct Vector_Traits: public Operator_Traits<ScalarType,SizeType,BackendType,Props...>
{};





template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUBackend,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::CPUBackend,Props...>
{
    using vector_type = std::vector<ScalarType>;
    using vector_type_reference = typename std::vector<ScalarType>::reference ;
    using scalar_type = typename std::vector<ScalarType>::value_type;
    using size_type = typename std::vector<ScalarType>::size_type;

};






#ifdef ChipSum_USE_KokkosKernels
template<typename ScalarType,typename SizeType,typename ...Props>
struct Vector_Traits<ScalarType,SizeType,
        ChipSum::Backend::KokkosKernels,Props...>:
        public Operator_Traits<ScalarType,SizeType,
        ChipSum::Backend::KokkosKernels,Props...>
{};
#endif





}
}
#endif // VECTOR_HPP
