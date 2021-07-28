/* * * * * * * * * * * * * * * * * * * * *
*   File:     backend.hpp
*   Author:   Li Kunyun
*   group:    CDCS-HPC
*   Time:     2021-07-28
* * * * * * * * * * * * * * * * * * * * * */

#ifndef __CHIPSUM_BACKEND_HPP__
#define __CHIPSUM_BACKEND_HPP__


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
