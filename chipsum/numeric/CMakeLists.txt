



list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/numeric_traits.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/sparse_matrix_types.h")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/sparse_matrix.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/dense_matrix.hpp")




list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl/vector_serial_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl/vector_kokkoskernels_impl.hpp")


list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl/crs_kokkoskernels_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl/crs_serial_impl.hpp")


list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl/densemat_blas_impl.hpp")

include_directories("../backend")

add_subdirectory("impl")

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
