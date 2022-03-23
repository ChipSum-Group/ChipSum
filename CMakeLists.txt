cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)

project(ChipSum VERSION 1.0)

configure_file(ChipSumConfig.h.in ChipSumConfig.h)

set(CMAKE_CXX_FLAGS "-Wall"${CMAKE_CXX_FLAGS})

set(CHIPSUM_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CHIPSUM_SOURCES "")

set(CHIPSUM_HEADERS "")


add_subdirectory("chipsum")


set(CMAKE_INCLUDE_DIRECTORIES_BEFORE on)

#option (ChipSum_ENABLE_KokkosKernels "Use provided Kokkos implementation" ON)

set(ChipSum_ENABLE_KokkosKernels ON)
if(ChipSum_ENABLE_KokkosKernels)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-build/lib)
        set(Kokkos_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-build/lib)
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-build/lib64)
        set(Kokkos_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-build/lib64)
    endif()

    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-kernels-build/lib)
        set(KokkosKernels_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-kernels-build/lib)
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-kernels-build/lib64)
        set(KokkosKernels_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpls/kokkos-kernels-build/lib64)
    endif()

    add_definitions("-DChipSum_USE_KokkosKernels=on")
    set(Kokkos_DIR ${Kokkos_LIB_DIR}/cmake/Kokkos)
    set(KokkosKernels_DIR ${KokkosKernels_LIB_DIR}/cmake/KokkosKernels)

    find_package(Kokkos REQUIRED)
    find_package(KokkosKernels REQUIRED)

    set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})

    list(APPEND CHIPSUM_LIBS  "Kokkos::kokkos")
    list(APPEND CHIPSUM_LIBS "Kokkos::kokkoskernels")

endif()

add_definitions("-fopenmp")

add_library(chipsum ${CHIPSUM_SOURCES} ${CHIPSUM_HEADERS})

list(APPEND CHIPSUM_LIBS chipsum)

install(TARGETS chipsum DESTINATION lib)
install(FILES ChipSum.hpp DESTINATION ${CMAKE_INSTALL_PREFIX})

set(CHIPSUM_ENABLE_TEST ON)
set(CHIPSUM_ENABLE_EXAMPLE OFF)
set(CHIPSUM_ENABLE_PERF_TEST ON)
set(CHIPSUM_ENABLE_UNIT_TEST ON)

if(CHIPSUM_ENABLE_TEST)
add_subdirectory("tests")
endif()
if(CHIPSUM_ENABLE_EXAMPLE)
add_subdirectory("examples")
endif()
if(CHIPSUM_ENABLE_PERF_TEST)
add_subdirectory("perf_test")
endif()
if(CHIPSUM_ENABLE_UNIT_TEST)
add_subdirectory("unit_test")
endif()

#"examples/chipsumAI/mnist.cpp" "examples/matrix/lu.cpp" "tests/test1.cpp" "perf_test/perf_test.cpp"
add_executable(ChipSum "examples/matrix/trtri.cpp")
install(TARGETS ChipSum  DESTINATION bin)
target_link_libraries(ChipSum PUBLIC ${CHIPSUM_LIBS})
target_include_directories(ChipSum PUBLIC ${PROJECT_BINARY_DIR})


