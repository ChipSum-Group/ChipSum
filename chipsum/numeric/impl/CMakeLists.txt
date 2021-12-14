# Core implementation code
add_subdirectory(cuda)
add_subdirectory(serial)
add_subdirectory(kokkoskernels)

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl_abstract.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector_impl.hpp")

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)

