cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(ChipSum VERSION 1.0)

configure_file(Testconfig.h.in Testconfig.h)


add_subdirectory("chipsum")
add_subdirectory("unit_test")

include_directories("chipsum")

add_executable(ChipSum
    "test.cpp"
    )

add_definitions("-DChipSum_USE_KokkosKernels")


target_link_libraries(ChipSum PUBLIC chipsum)

target_include_directories(ChipSum PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

