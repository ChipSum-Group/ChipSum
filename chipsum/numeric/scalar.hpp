/*
 * @Author: your name
 * @Date: 2021-08-09 12:27:29
 * @LastEditTime: 2021-08-09 16:31:33
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: scalar.hpp
 */
#ifndef __CHIPSUM_NUMERIC_SCALAR_HPP_
#define __CHIPSUM_NUMERIC_SCALAR_HPP_

#include "../chipsum_macro.h"
#include "numeric_traits.hpp"
#include "impl/scalar_kokkoskernels_impl.hpp"
#include "impl/scalar_serial_impl.hpp"
#include "../backend/backend.hpp"
#include "vector.hpp"

namespace ChipSum
{
    namespace Numeric
    {

        template <typename... Props>
        class Scalar;

        template <typename ScalarType, typename SizeType, typename BackendType, typename... Props>
        class Scalar<ScalarType, SizeType, BackendType, Props...>
        {

        public:
            using traits = Scalar_Traits<ScalarType, SizeType, BackendType, Props...>;

            using scalar_type = typename traits::scalar_type;

            using scalar_type_reference = typename std::add_lvalue_reference<scalar_type>::type;
            using const_scalar_type_reference = typename std::add_const<scalar_type_reference>::type;

            using device_scalar_type = typename traits::device_scalar_value_type;
            using device_scalar_type_reference = typename std::add_lvalue_reference<device_scalar_type>::type;
            ;

        private:
            scalar_type __data;

        public:
            CHIPSUM_DECLARED_FUNCTION Scalar()
            {
                ChipSum::Numeric::Impl::Scalar::Create<ScalarType, SizeType>(__data);
            }

            CHIPSUM_DECLARED_FUNCTION Scalar(const ScalarType s)
            {
                ChipSum::Numeric::Impl::Scalar::Create<ScalarType, SizeType>(s, __data);
            }


            CHIPSUM_FUNCTION_INLINE const_scalar_type_reference GetData(){return __data;}

            CHIPSUM_FUNCTION_INLINE Scalar operator=(ScalarType &&s)
            {
                ChipSum::Numeric::Impl::Scalar::DeepCopy<ScalarType, SizeType>(s, __data);
            }

            CHIPSUM_FUNCTION_INLINE const ScalarType operator()(){
                return ChipSum::Numeric::Impl::Scalar::GetItem<ScalarType,SizeType>(__data);
            }
        };


        template<typename ScalarType,typename SizeType,typename BackendType,typename ...Props>
        CHIPSUM_FUNCTION_INLINE Vector<ScalarType,SizeType,BackendType> 
        operator*(Scalar<ScalarType,SizeType,BackendType>& s,Vector<ScalarType,SizeType,BackendType>& v)
        {
            Vector<ScalarType,SizeType,BackendType> r(v.GetSize());
            ChipSum::Numeric::Impl::Scalar::Mult<ScalarType,SizeType>(s.GetData(),v.GetData(),r.GetData());
            return r;
        }



    } // End namespace Numeric
} // End namespace ChipSum

typedef ChipSum::Numeric::Scalar<double, std::size_t, ChipSum::Backend::DefaultBackend> Scalar;

#endif