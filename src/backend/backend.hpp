#ifndef BACKEND_HPP
#define BACKEND_HPP


namespace ChipSum{
namespace Backend{

template<typename ...Props>
struct Backend_Traits;

template<typename ...Props>
struct Builtin:public Backend_Traits<Props...>{/*TODO*/};


template<typename ...Props>
struct KokkosKernels:public Backend_Traits<Props...>{/*TODO*/};

}
};


#endif // BACKEND_HPP
