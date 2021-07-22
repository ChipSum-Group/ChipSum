cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

project(ChipSum VERSION 1.0)

configure_file(ChipSumConfig.h.in ChipSumConfig.h)


set(CHIPSUM_SOURCES "")

set(CHIPSUM_HEADERS "")

add_subdirectory("unit_test")
add_subdirectory("chipsum")


#set(ChipSum_USE_KokkosKernels "yes")
#set(ChipSum_Kokkos_DIR "/home/lky/code/git/kokkos-build")
#set(ChipSum_KokkosKernels_DIR "/home/lky/code/git/kokkos-kernels-build")


if(ChipSum_USE_KokkosKernels)
    set(KokkosKernels_DIR ${ChipSum_KokkosKernels_DIR}/lib/cmake/KokkosKernels)
    set(Kokkos_DIR ${ChipSum_Kokkos_DIR}/lib/cmake/Kokkos)

    find_package(Kokkos REQUIRED)
    find_package(KokkosKernels REQUIRED)
    set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})
    list(APPEND CHIPSUM_LIBS  "Kokkos::kokkos")
    list(APPEND CHIPSUM_LIBS "Kokkos::kokkoskernels")

endif()












add_definitions("-fopenmp")

include_directories("${ChipSum_Kokkos_DIR}/include")
include_directories("${ChipSum_KokkosKernels_DIR}/include")

add_library(chipsum ${CHIPSUM_SOURCES} ${CHIPSUM_HEADERS})

list(APPEND CHIPSUM_LIBS chipsum)

install(TARGETS chipsum DESTINATION lib)
install(FILES ${CHIPSUM_HEADERS}  DESTINATION include)


add_executable(ChipSum
    "test.cpp"
    )


target_link_libraries(ChipSum PUBLIC ${CHIPSUM_LIBS}  )

target_include_directories(ChipSum PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )



