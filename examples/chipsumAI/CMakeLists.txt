cmake_minimum_required(VERSION 3.16)

project(untitled LANGUAGES CXX)



set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)

add_definitions("-DChipSum_USE_KokkosKernels=On")

include_directories(${ChipSum_DIR}/include)
link_directories(${ChipSum_DIR}/lib)

set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})

add_executable(mnist mnist.cpp)

target_link_libraries(mnist Kokkos::kokkos Kokkos::kokkoskernels chipsum)

