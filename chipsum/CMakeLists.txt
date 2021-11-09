





set(SOURCES "")
set(HEADERS "")

add_subdirectory(common)
add_subdirectory(numeric)
add_subdirectory(backend)
#add_subdirectory(distribution)



include_directories(".")

set(CHIPSUM_SOURCES ${CHIPSUM_SOURCES} ${SOURCES} PARENT_SCOPE)
set(CHIPSUM_HEADERS ${CHIPSUM_HEADERS} ${HEADERS} PARENT_SCOPE)



#message(${Kokkos_CXX_COMPILER})



