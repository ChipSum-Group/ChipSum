# Core implementation code



list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/scalar_kokkoskernels_impl.hpp")


list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_dot_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_scal_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_nrm1_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_nrm2_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_nrminf_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_kokkoskernels_axpby_impl.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/csr_kokkoskernels_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/csr_kokkoskernels_spmv_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/csr_kokkoskernels_spgemm_impl.hpp")


list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_scal_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_gemv_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_gemm_impl.hpp")

# list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_relu_impl.hpp")
# list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_activation_impl.hpp")
# list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_norm_impl.hpp")
# list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_kokkoskernels_dense_impl.hpp")

install(FILES ${HDRS} DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum/numeric/impl/kokkoskernels)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
