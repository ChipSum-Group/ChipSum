#ifndef __CHIPSUM_CSR_SERIAL_SPMV_IMPL_HPP__
#define __CHIPSUM_CSR_SERIAL_SPMV_IMPL_HPP__

#include "../../../chipsum_macro.h"

namespace  ChipSum{
namespace  Numeric{
namespace Impl {
namespace Sparse {
template <typename ValueType,typename ST>
// b=Ax 一种常用的SpMV接口
CHIPSUM_FUNCTION_INLINE void
spmv(ST &A,::std::vector<ValueType> &x, ::std::vector<ValueType> &b) {


  ::std::size_t M = A.graph.row_map.size()-1;

  for (::std::size_t i = 0; i < M; ++i)
    b[i] = 0;

  for (::std::size_t i = 0; i < M; ++i) {
    ::std::size_t start = A.graph.row_map[i];
    ::std::size_t end = A.graph.row_map[i + 1];
    for (::std::size_t j = 0; j < end - start; ++j) {
      b[i] += A.vals[start + j] * x[A.graph.col_map[start + j]];
    }
  }
}

template <typename ValueType,
          typename ST,
          typename AlphaT,
          typename BetaT
          >
// b = beta*b+alpha*A*x 完整的SpMV
CHIPSUM_FUNCTION_INLINE void
spmv( ST &A,
     ::std::vector<ValueType> &x,
      ::std::vector<ValueType> &b,
      const AlphaT& alpha,
      const BetaT& beta
      ) {

  for (::std::size_t i = 0; i < b.size(); ++i) {
    ::std::size_t start = A.graph.row_map[i];
    ::std::size_t end = A.graph.row_map[i + 1];
    for (::std::size_t j = 0; j < end - start; ++j) {
      b[i] += beta * b[i] +
              alpha * A.vals[start + j] * x[A.graph.col_map[start + j]];
    }
  }
}


}
}
}
}
#endif // CSR_SERIAL_SPMV_IMPL_HPP
