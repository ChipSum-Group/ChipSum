#ifndef __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__
#define __CHIPSUM_NUMERIC_SPARSE_MATRIX_HPP__

#include <vector>
#include <type_traits>

//#include <iostream>
//using namespace std;

#include "numeric_traits.hpp"
#include "impl/vector_serial_impl.hpp"
#include "impl/vector_kokkoskernels_impl.hpp"




namespace ChipSum {
namespace Numeric {

template<typename ...Props>
class SparseMatrix;

template<typename ScalarType,typename SizeType,typename SpFormat,typename BackendType,typename ...Props>
class SparseMatrix<ScalarType,SizeType,SpFormat,BackendType,Props...>{

public:







};






} // End namespace Numeric
} // End namespace ChipSum



#endif // SPARSEMATRIX_HPP
