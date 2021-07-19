#ifndef BACKEND_HPP
#define BACKEND_HPP


namespace ChipSum{
namespace Backend{



struct Backend{};


struct Builtin:public Backend{/*TODO*/};


typedef  Builtin CPUBackend;






struct KokkosKernels:public Backend{/*TODO*/};

#ifdef ChipSum_USE_KokkosKernels
typedef KokkosKernels DefaultBackend;

#else
typedef Builtin DefaultBackend;

#endif
}
};


#endif // BACKEND_HPP
