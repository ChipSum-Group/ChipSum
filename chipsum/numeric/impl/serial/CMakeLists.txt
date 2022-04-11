# Core implementation code

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/scalar_serial_impl.hpp")


list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_scal_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_dot_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_axpby_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_nrm1_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_nrm2_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_serial_nrminf_impl.hpp")




list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_serial_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_serial_gemm_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_serial_gemv_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/densemat_serial_scal_impl.hpp")





list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/csr_serial_impl.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/csr_serial_spmv_impl.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/coo_serial_impl.hpp")


install(FILES ${HDRS} DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum/numeric/impl/serial)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
