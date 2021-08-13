/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-13 10:21:45
 */


#ifndef __CHIPSUM_DENSE_MATRIX_HPP__
#define __CHIPSUM_DENSE_MATRIX_HPP__

#include "impl/densemat_kokkoskernels_impl.hpp"
#include "impl/densemat_serial_impl.hpp"
#include "numeric_traits.hpp"
#include "scalar.hpp"
#include "vector.hpp"

namespace ChipSum {
namespace Numeric {

template <typename... Props> class DenseMatrix;

template <typename ScalarType, typename SizeType, typename BackendType,
          typename... Props>
class DenseMatrix<ScalarType, SizeType, BackendType, Props...> {

public:
  using traits =
      DenseMatrix_Traits<ScalarType, SizeType, BackendType, Props...>;
  using matrix_type = typename traits::matrix_type;
  using matrix_type_reference =
      typename std::add_lvalue_reference<matrix_type>::type;
  using const_matrix_type_reference =
      typename std::add_const<matrix_type_reference>::type;

  using size_type = typename traits::size_type;
  using const_size_type = typename traits::const_size_type;

  using vector_type =
      ChipSum::Numeric::Vector<ScalarType, SizeType, BackendType, Props...>;

private:
  matrix_type __data;
  size_type __nrow;
  size_type __ncol;

public:
  /**
   * @description:
   * @param {size_type} M
   * @param {size_type} N
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION DenseMatrix(const_size_type M, const_size_type N)
      : __nrow(M), __ncol(N) {
    ChipSum::Numeric::Impl::DenseMat::Create<ScalarType, SizeType>(M, N,
                                                                   __data);
  }

  /**
   * @description:
   * @param {size_type} M
   * @param {size_type} N
   * @param {ScalarType*} src
   * @return {*}
   */
  CHIPSUM_DECLARED_FUNCTION DenseMatrix(const_size_type M, const_size_type N,
                                        ScalarType *src)
      : __nrow(M), __ncol(N) {
    ChipSum::Numeric::Impl::DenseMat::Create<ScalarType, SizeType>(M, N,
                                                                   __data);
    ChipSum::Numeric::Impl::DenseMat::Fill<ScalarType, SizeType>(M, N, src,
                                                                 __data);
                                                                 Kokkos::fence();
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE const_matrix_type_reference GetData() {
    return __data;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE size_type GetRowNum() { return __nrow; }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE size_type GetColNum() { return __ncol; }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE DenseMatrix operator*(DenseMatrix &m) {
    DenseMatrix ret(__nrow, m.GetColNum());
    ChipSum::Numeric::Impl::DenseMat::Mult<ScalarType, SizeType>(
        __nrow,m.GetColNum(),m.GetRowNum(), __data,m.GetData(), ret.GetData());
    return ret;
  }

  // /**
  //  * @description:
  //  * @param {*}
  //  * @return {*}
  //  */
  // CHIPSUM_FUNCTION_INLINE DenseMatrix operator*=(DenseMatrix &m) {
    
  //   ChipSum::Numeric::Impl::DenseMat::Mult<ScalarType, SizeType>(
  //       __nrow,m.GetColNum(),m.GetRowNum(), __data,m.GetData(), ret.GetData());
  //   return ret;
  // }

  template <typename... Args>
  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE vector_type operator*(vector_type &v) {
    vector_type ret(__ncol);
    ChipSum::Numeric::Impl::DenseMat::Mult<ScalarType, SizeType>(
        __nrow, __ncol, __data, v.GetData(), ret.GetData());
    return ret;
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE matrix_type operator*=(ScalarType s) {
    ChipSum::Numeric::Impl::DenseMat::Scal<ScalarType, SizeType>(s, __data);
    return *this;
  }
  /**
   * @description: 
   * @param {*}
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE ScalarType &operator()(const_size_type i,
                                                 const_size_type j) {
    return ChipSum::Numeric::Impl::DenseMat::GetItem<ScalarType, SizeType>(
        i, j, __nrow, __ncol, __data);
  }

    /**
   * @description: 
   * @param {*}
   * @return {*}
   * @author: Li Kunyun
   */
  CHIPSUM_FUNCTION_INLINE const ScalarType &operator()(const_size_type i,
                                                 const_size_type j) const{
    return ChipSum::Numeric::Impl::DenseMat::GetItem<ScalarType, SizeType>(
        i, j, __nrow, __ncol, __data);
  }

  /**
   * @description:
   * @param {*}
   * @return {*}
   */
  CHIPSUM_FUNCTION_INLINE void Print(std::ostream &out = std::cout) {
    ChipSum::Numeric::Impl::DenseMat::Print<ScalarType, SizeType>(
        __nrow, __ncol, __data, out);
  }
};

} // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::DenseMatrix<double, std::size_t,
                                      ChipSum::Backend::DefaultBackend>
    Matrix;

#endif // __CHIPSUM_DENSE_MATRIX_HPP__
