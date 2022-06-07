





set(SOURCES "")
set(HEADERS "")



add_subdirectory(common)
add_subdirectory(numeric)
add_subdirectory(backend)
add_subdirectory(solver)
#add_subdirectory(distribution)



include_directories(".")

set(CHIPSUM_SOURCES ${CHIPSUM_SOURCES} ${SOURCES} PARENT_SCOPE)
set(CHIPSUM_HEADERS ${CHIPSUM_HEADERS} ${HEADERS} PARENT_SCOPE)

install(FILES ChipSumConfig.h DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum)
install(FILES chipsum_macro.h DESTINATION  ${CMAKE_INSTALL_PREFIX}/include/chipsum)


#message(${Kokkos_CXX_COMPILER})



