#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__

#include <vector>
#include <type_traits>

//#include <iostream>
//using namespace std;

#include "numeric_traits.hpp"
#include "impl/crs_kokkoskernels_impl.hpp"





namespace ChipSum {
namespace Numeric {

template<typename ...Props>
class SparseMatrix;

template<typename ScalarType,typename SizeType,typename SpFormat,typename BackendType,typename ...Props>
class SparseMatrix<ScalarType,SizeType,SpFormat,BackendType,Props...>{

public:
    using traits = Sparse_Traits<ScalarType,SizeType,BackendType,Props...>;


    using vector_type = typename traits::vector_type;
    using size_type = typename traits::size_type;
    using size_type_reference = typename std::add_lvalue_reference<size_type>::type;
    using const_size_type_reference = typename std::add_const<size_type_reference>::type;

    using vector_type_reference = typename std::add_lvalue_reference<vector_type>::type;
    using const_vector_type_reference = typename std::add_const<vector_type_reference>::type;






};






} // End namespace Numeric
} // End namespace ChipSum



#endif // SPARSEMATRIX_HPP
