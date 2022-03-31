cmake_minimum_required(VERSION 3.16)

project(untitled LANGUAGES CXX)



set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(EXISTS ${ChipSum_DIR}/tpls/kokkos-build/lib)
    set(Kokkos_LIB_DIR ${ChipSum_DIR}/tpls/kokkos-build/lib)
elseif(EXISTS ${ChipSum_DIR}/tpls/kokkos-build/lib64)
    set(Kokkos_LIB_DIR ${ChipSum_DIR}/tpls/kokkos-build/lib64)
endif()

if(EXISTS ${ChipSum_DIR}/tpls/kokkos-kernels-build/lib)
    set(KokkosKernels_LIB_DIR ${ChipSum_DIR}/tpls/kokkos-kernels-build/lib)
elseif(EXISTS ${ChipSum_DIR}/tpls/kokkos-kernels-build/lib64)
    set(KokkosKernels_LIB_DIR ${ChipSum_DIR}/tpls/kokkos-kernels-build/lib64)
endif()

set(Kokkos_DIR ${Kokkos_LIB_DIR}/cmake/Kokkos)
set(KokkosKernels_DIR ${KokkosKernels_LIB_DIR}/cmake/KokkosKernels)

find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)

add_definitions("-DChipSum_USE_KokkosKernels=On")

include_directories(${ChipSum_DIR}/include)
link_directories(${ChipSum_DIR}/lib)

set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})

add_executable(mnist mnist.cpp)

target_link_libraries(mnist Kokkos::kokkos Kokkos::kokkoskernels chipsum)

