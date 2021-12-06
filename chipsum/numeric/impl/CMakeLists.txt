# Core implementation code
add_subdirectory(cuda)
add_subdirectory(serial)
add_subdirectory(kokkoskernels)

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/impl_abstract.hpp")

#install(FILES ${HDRS} DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum/numeric/impl)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)

