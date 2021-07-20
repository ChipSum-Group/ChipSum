cmake_minimum_required(VERSION 3.5)

set(SOURCES "")
set(HEADERS "")

add_subdirectory(common)
add_subdirectory(numeric)
add_subdirectory(backend)



add_library(chipsum ${SOURCES} ${HEADERS})


include_directories(
    common
    numeric
    backend
    )

install(TARGETS chipsum DESTINATION lib)
install(FILES ${HEADERS}  DESTINATION include)
