#ifndef BACKEND_HPP
#define BACKEND_HPP


namespace ChipSum{
namespace Backend{



struct Backend{};


struct BuiltinSerial:public Backend{/*TODO*/};


typedef  BuiltinSerial CPUSerialBackend;






struct KokkosKernels:public Backend{/*TODO*/};

#ifdef ChipSum_USE_KokkosKernels
typedef KokkosKernels DefaultBackend;

#else
typedef BuiltinSerial DefaultBackend;

#endif
}
};


#endif // BACKEND_HPP
