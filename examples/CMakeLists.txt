# Minimal required: cmake 3.15
add_subdirectory("vector")
add_subdirectory("sparse")
add_subdirectory("matrix")
add_subdirectory("tensor")
add_subdirectory("chipsumAI")
file(GLOB testfiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(mainfile IN LISTS testfiles)
    #get file name without directory
    get_filename_component(mainname ${mainfile} NAME_WE)
    add_executable(${mainname} ${mainfile})
    target_link_libraries(${mainname} ${CHIPSUM_LIBS})
    target_include_directories(${mainname} PUBLIC ${PROJECT_BINARY_DIR})
endforeach()
