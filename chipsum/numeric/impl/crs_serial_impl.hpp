/* * * * * * * * * * * * * * * * * * * * *
*   File:     crs_kokkoskernels_impl.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_CRS_SERIAL_IMPL_HPP__
#define __CHIPSUM_CRS_SERIAL_IMPL_HPP__


#include <vector>

#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"
#include "../../chipsum_macro.h"


namespace ChipSum{
namespace Numeric {

template<typename SizeType,typename ...Props>
struct StaticGraph{
    std::vector<SizeType> row_map;
    std::vector<SizeType> col_map;
};


template<typename ScalarType,typename SizeType,typename ...Props>
struct CrsFormat{
    std::vector<ScalarType> vals;
    StaticGraph<SizeType> graph;

};


template<typename ScalarType,typename SizeType,typename ...Props>
struct Sparse_Traits<ScalarType,SizeType,SparseTypes::Csr,ChipSum::Backend::BuiltinSerial,Props...>
        : public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::BuiltinSerial,Props...>{


    using matrix_format_type = CrsFormat<ScalarType,SizeType>;

    using graph_type = StaticGraph<SizeType>;
    using row_map_type = std::vector<SizeType>;
    using col_map_type = std::vector<SizeType>;
    using matrix_values_type = std::vector<ScalarType>;


};




namespace Impl {

namespace Sparse {

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Fill：//TODO
 * @param nrows
 * @param ncols
 * @param annz
 * @param row_map
 * @param col_map
 * @param values
 * @param A
 */
CHIPSUM_FUNCTION_INLINE void Create(CrsFormat<ScalarType,SizeType>& A,
                                    const SizeType nrows,
                                    const SizeType ncols,
                                    const SizeType annz,
                                    SizeType* row_map,
                                    SizeType* col_map,
                                    ScalarType* values

                                    )
{
    A.vals = std::vector<ScalarType>(values,values+annz);
    A.graph.row_map = std::vector<SizeType>(row_map,row_map+nrows+1);
    A.graph.col_map = std::vector<SizeType>(col_map,col_map+annz);

}

template <typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Create
 * @param A
 * @param row_map_size
 * @param col_map_size
 */
CHIPSUM_FUNCTION_INLINE void Create(CrsFormat<ScalarType,SizeType>& A,
                                    const std::size_t row_map_size,
                                    const std::size_t col_map_size
                                    )
{
    A.vals.resize(col_map_size);
    A.graph.row_map.resize(col_map_size);
    A.graph.col_map.resize(col_map_size);
}

template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Mult
 * @param A
 * @param x
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void Mult(CrsFormat<ScalarType,SizeType>& A,
                                  std::vector<ScalarType>& x,
                                  std::vector<ScalarType>& b)
{

    for(std::size_t i=0;i<b.size();++i) b[i]=0;

    for(std::size_t i=0;i<b.size();++i)
    {
        std::size_t start = A.graph.row_map[i];
        std::size_t end   = A.graph.row_map[i+1];
        for(std::size_t j=0;j<end-start;++j)
        {
            b[i] += A.vals[start+j]*x[A.graph.col_map[start+j]];


        }

    }
}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Mult
 * @param alpha
 * @param A
 * @param x
 * @param beta
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void Mult(
        ScalarType alpha,
        CrsFormat<ScalarType,SizeType>& A,
        std::vector<ScalarType>& x,
        ScalarType beta,
        std::vector<ScalarType>& b)
{


    for(std::size_t i=0;i<b.size();++i)
    {
        std::size_t start = A.graph.row_map[i];
        std::size_t end   = A.graph.row_map[i+1];
        for(std::size_t j=0;j<end-start;++j)
        {
            b[i] += beta*b[i] + alpha*A.vals[start+j]*x[A.graph.col_map[start+j]];
        }
    }
}

template <typename ScalarType,typename SizeType,typename ...Props>
CHIPSUM_FUNCTION_INLINE void Gauss_seidel(/*...*/){
    //TODO
}


} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
