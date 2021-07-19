cmake_minimum_required(VERSION 3.5)

project(numeric)

list(APPEND SRCS "${PROJECT_SOURCE_DIR}/operator.cpp")
list(APPEND HDRS "${PROJECT_SOURCE_DIR}/operator.hpp")


list(APPEND SRCS "${PROJECT_SOURCE_DIR}/vector.cpp")
list(APPEND HDRS "${PROJECT_SOURCE_DIR}/vector.hpp")


list(APPEND SRCS "${PROJECT_SOURCE_DIR}/impl/vector_impl.cpp")
list(APPEND HDRS "${PROJECT_SOURCE_DIR}/impl/vector_impl.hpp")


include_directories("../backend")


set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
