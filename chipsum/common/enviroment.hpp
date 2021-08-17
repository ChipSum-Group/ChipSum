/*
 * @Description:
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-17 09:30:04
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-17 15:56:04
 */

#ifndef __CHIPSUM_ENVIROMENT_HPP__
#define __CHIPSUM_ENVIROMENT_HPP__

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include <Kokkos_Core.hpp>
#endif
namespace ChipSum {
namespace Common {

void Init(int &argc, char *argv[]) {
#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
  Kokkos::initialize(argc, argv);
#endif
}

void Finalize() {
#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
  Kokkos::finalize();
#endif
}

} // namespace Common
} // namespace ChipSum

#endif
