///
/// \file     scalar.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-01
/// \brief    标量用户接口，主要是为了衔接类似点积一类
///           操作的Device端实现。
///


#ifndef __CHIPSUM_NUMERIC_SCALAR_HPP__
#define __CHIPSUM_NUMERIC_SCALAR_HPP__

#include "../backend/backend.hpp"
#include "../chipsum_macro.h"

#include "numeric_traits.hpp"

#include "impl/serial/scalar_serial_impl.hpp"

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include "impl/kokkoskernels/scalar_kokkoskernels_impl.hpp"
// #include "impl/scalar_kokkos_impl.hpp"
#endif



namespace ChipSum {
namespace Numeric {

template <typename... Props>
class Scalar {

public:
    using traits = Scalar_Traits<Props...>;

    using scalar_type = typename traits::scalar_type;
    using scalar_type_ref =
    typename ::std::add_lvalue_reference<scalar_type>::type;
    using const_scalar_type_ref =
    typename ::std::add_const<scalar_type_ref>::type;


    using host_scalar_type = typename traits::host_type;
    using host_scalar_type_ref =
    typename std::add_lvalue_reference<host_scalar_type>::type;
    using const_host_scalar_type_ref =
    typename ::std::add_const<host_scalar_type_ref>::type;


    ///
    /// \brief like double or something
    ///
    using value_type = typename traits::value_type;
    using value_type_ref =
    typename std::add_lvalue_reference<value_type>::type;
    using const_value_type_ref =
    typename ::std::add_const<value_type_ref>::type;

    private:
    scalar_type __data;

public:



    CHIPSUM_DECLARED_FUNCTION
    ///
    /// \brief Scalar 构造函数
    ///
    Scalar() {
        ChipSum::Numeric::Impl::Scalar::create(__data);
    }


    CHIPSUM_DECLARED_FUNCTION
    ///
    /// \brief Scalar 构造函数
    ///
    Scalar(const value_type& s) {
        ChipSum::Numeric::Impl::Scalar::create(__data, s);
    }


    CHIPSUM_FUNCTION_INLINE
    ///
    /// \brief GetData 获取后端底层，如Kokkos::View<double>
    /// \return 后端数据引用
    ///
    const_scalar_type_ref GetData() {
        return __data;
    }


    CHIPSUM_FUNCTION_INLINE
    typename ::std::add_lvalue_reference<Scalar>::type
    operator=(const value_type& s) {
        ChipSum::Numeric::Impl::Scalar::deep_copy( __data,s);
        return *this;
    }

    CHIPSUM_FUNCTION_INLINE
    ///
    /// \brief operator ()
    /// \return
    ///
    value_type operator()() {
        value_type r;
        ChipSum::Numeric::Impl::Scalar::get_item(
                    __data,r);
        return r;
    }

    CHIPSUM_FUNCTION_INLINE
    ///
    /// \brief operator ScalarType
    ///
    operator value_type()  {
        value_type r;
        ChipSum::Numeric::Impl::Scalar::get_item(
                    __data,r);
        return r;
    }

    CHIPSUM_FUNCTION_INLINE
    ///
    /// \brief Print
    /// \param out
    ///
    void Print(std::ostream &out = std::cout) {
        ChipSum::Numeric::Impl::Scalar::print(__data, out);
    }
};

} // End namespace Numeric
} // End namespace ChipSum




typedef ChipSum::Numeric::Scalar<CSFloat,ChipSum::Backend::DefaultBackend>
Scalar;

#endif
