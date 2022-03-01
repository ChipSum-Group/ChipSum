



# list(APPEND SRCS "${CMAKE_CURRENT_SOURCE_DIR}/enviroment.cpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/enviroment.hpp")



list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/bmp_writer.h")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/bmp_writer.cpp")

list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/png_writer.hpp")
list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/coo_reader.hpp")
# list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/bmp_writer.cpp")


list(APPEND HDRS "${CMAKE_CURRENT_SOURCE_DIR}/data_types.hpp")

install(FILES ${HDRS} DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum/common)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
