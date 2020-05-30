cmake_minimum_required(VERSION 3.16)
project(LocallyInjectiveMappings)

set(CMAKE_CXX_STANDARD 11)

set(LIB_SURFACE_MESH /usr/local/lib/libSurfaceMesh.1.0.dylib)

include_directories(/usr/local/include)

find_package(Eigen3 3.3.0)
if (Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(can't find eigen3)
    return()
endif()

aux_source_directory(. DIR_SRC)

add_executable(${PROJECT_NAME} ${DIR_SRC})

target_link_libraries(${PROJECT_NAME} ${LIB_SURFACE_MESH})