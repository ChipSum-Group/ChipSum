#ifndef __VECTOR_IMPL_HPP__
#define __VECTOR_IMPL_HPP__


#include "../vector.hpp"
namespace ChipSum{
namespace Numeric {

namespace Impl {



template<typename ...Props>
void Dot(Props...);

template <typename ScalarType,
          typename SizeType,
          typename Backend,
          typename ...Props
          >
void Dot(ChipSum::Numeric::Vector_Traits<ScalarType,SizeType,Backend,Props...> a){
    cout<<"aaa"<<endl;
}


}

}

}

#endif
