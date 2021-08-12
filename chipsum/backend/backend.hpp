/*
 * @Description: 后端模板参数
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 12:20:42
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 15:22:01
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

// 适合用于小矩阵运算，考虑后面为它设计默认编译选项，
// 形成BlasMatrix。
struct OpenBlas:public BackendBase{/*TODO*/};

#ifdef ChipSum_USE_KokkosKernels || ChipSum_USE_KokkosKernels64
typedef KokkosKernels DefaultBackend;

#else
typedef KokkosKernels DefaultBackend;

#endif
}
}


#endif // BACKEND_HPP
