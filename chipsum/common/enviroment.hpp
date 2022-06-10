/*
 * @Description:
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-17 09:30:04
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2022-06-09 16:38:59
 */

#ifndef __CHIPSUM_ENVIROMENT_HPP__
#define __CHIPSUM_ENVIROMENT_HPP__

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#include <Kokkos_Core.hpp>
#endif
namespace ChipSum {
namespace Common {


///
/// \brief Init ChipSum环境初始化
/// \param argc 命令行参数argc
/// \param argv 命令行参数argv
///
void Init(int &argc, char *argv[]) {
#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
  Kokkos::initialize(argc, argv);
#endif
}


///
/// \brief Finalize 安全退出ChipSum环境
///
void Finalize() {
#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
  Kokkos::finalize();
#endif
}

} // namespace Common
} // namespace ChipSum

#endif
