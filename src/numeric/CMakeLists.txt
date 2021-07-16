cmake_minimum_required(VERSION 3.5)

project(numeric)

list(APPEND SRCS "${PROJECT_SOURCE_DIR}/operator.cpp")
list(APPEND HDRS "${PROJECT_SOURCE_DIR}/operator.hpp")


include_directories("../backend")


set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
