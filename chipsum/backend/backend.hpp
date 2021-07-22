#ifndef BACKEND_HPP
#define BACKEND_HPP


namespace ChipSum{
namespace Backend{



struct BackendBase{};


struct BuiltinSerial:public BackendBase{/*TODO*/};


typedef  BuiltinSerial CPUSerialBackend;



struct KokkosKernels:public BackendBase{/*TODO*/};

#ifdef ChipSum_USE_KokkosKernels
typedef KokkosKernels DefaultBackend;

#else
typedef BuiltinSerial DefaultBackend;

#endif
}
}


#endif // BACKEND_HPP
