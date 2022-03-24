



list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/numeric_traits.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/scalar.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/vector.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/sparse_matrix_types.h")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/sparse_matrix.hpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/dense_matrix.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/tensor.hpp")




add_subdirectory("impl")

install(FILES ${HDRS} DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum/numeric)



set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
