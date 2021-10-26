/*
 * @Description: 后端模板参数
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-16 10:30:02
 */


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
typedef KokkosKernels DefaultBackend;

#else

typedef Serial DefaultBackend;

#endif
}
}


#endif // BACKEND_HPP
