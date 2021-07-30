/* * * * * * * * * * * * * * * * * * * * *
*   File:     crs_kokkoskernels_impl.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
#define __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__


#include <vector>

#include "../numeric_traits.hpp"
#include "../sparse_matrix_types.h"
#include "../../chipsum_macro.h"


static int spm_name = 0;

namespace ChipSum{
namespace Numeric {

template<typename ScalarType,typename ...Props>
struct StaticGraph{
    std::vector<typename std::vector<ScalarType>::size_type> row_map;
    std::vector<typename std::vector<ScalarType>::size_type> col_map;
};


template<typename ScalarType,typename ...Props>
struct CrsFormat{
    std::vector<ScalarType> vals;
    StaticGraph<ScalarType> graph;

};


template<typename ScalarType,typename SizeType,typename ...Props>
struct Sparse_Traits<ScalarType,SizeType,Csr,ChipSum::Backend::BuiltinSerial,Props...>
        : public Operator_Traits<ScalarType,SizeType,ChipSum::Backend::BuiltinSerial,Props...>{


    using matrix_format_type = CrsFormat<ScalarType>;

    using graph_type = StaticGraph<ScalarType>;
    using row_map_type = std::vector<typename std::vector<ScalarType>::size_type>;
    using col_map_type = std::vector<typename std::vector<ScalarType>::size_type>;
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
CHIPSUM_FUNCTION_INLINE void Fill(
        CrsFormat<ScalarType>& A,
        const SizeType nrows,
        const SizeType ncols,
        const SizeType annz,
        SizeType* row_map,
        SizeType* col_map,
        ScalarType* values

        )
{
    A.vals = std::vector<ScalarType>(values,values+annz);
    A.graph.row_map = std::vector<decltype (A.graph.row_map)>(row_map,row_map+nrows+1);
    A.graph.col_map = std::vector<decltype (A.graph.col_map)>(col_map,col_map+annz);

}


template<typename ScalarType,typename SizeType,typename ...Props>
/**
 * @brief Mult
 * @param A
 * @param x
 * @param b
 */
CHIPSUM_FUNCTION_INLINE void Mult(
        CrsFormat<ScalarType>& A,
        std::vector<ScalarType>& x,
        std::vector<ScalarType>& b)
{

    for(int i=0;i<b.size();++i) b[i]=0.0;

    for(size_t i=0;i<b.size();++i)
    {
        size_t start = A.graph.row_map[i];
        size_t end   = A.graph.row_map[i+1];
        for(size_t j=start,j<end;++j)
        {
            b[i] += A.vals[j]*x[j];
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
        CrsFormat<ScalarType>& A,
        std::vector<ScalarType>& x,
        ScalarType beta,
        std::vector<ScalarType>& b)
{
    //TODO
}


} // End namespace Sparse
} // End namespace Impl
} // End namespace Numeric
} // End namespace ChipSum

#endif // __CHIPSUM_CRS_KOKKOSKERNELS_IMPL_HPP__
