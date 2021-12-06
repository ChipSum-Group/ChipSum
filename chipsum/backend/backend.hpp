///
/// \file     backend.hpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-12-02
/// \brief    %stuff%
///

#ifndef __CHIPSUM_BACKEND_HPP__
#define __CHIPSUM_BACKEND_HPP__


namespace ChipSum{
namespace Backend{

struct BackendBase{};

struct Serial:public BackendBase{/*TODO*/};

typedef  Serial CPUSerial;

struct KokkosKernels:public BackendBase{/*TODO*/};

struct Kokkos:public BackendBase{/*TODO*/};

struct Cuda:public BackendBase{/*TODO*/};

// 适合用于小矩阵运算，考虑后面为它设计默认编译选项，
// 形成BlasMatrix。
struct OpenBlas:public BackendBase{/*TODO*/};

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
typedef Serial DefaultBackend;

#else

typedef Serial DefaultBackend;

#endif
}
}


#endif // BACKEND_HPP
