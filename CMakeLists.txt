cmake_minimum_required(VERSION 3.0)

project(mesh-analyzer)

set(CMAKE_CXX_STANDARD_REQUIRED 11)
set(CMAKE_CXX_FLAGS '-std=c++11')
set(CMAKE_INSTALL_PREFIX /usr/local/bin)

# Next lines needed for building all Qt projects
find_package( Qt4 REQUIRED COMPONENTS QtXml)
include( ${QT_USE_FILE} )
add_definitions( ${QT_DEFINITIONS} )
include_directories( ${CMAKE_BINARY_DIR} )

find_package(OpenGL REQUIRED COMPONENTS GLU)
include_directories("${OPENGL_INCLUDE_DIR}")

include_directories(.)

message("${QT_LIBRARIES}")

add_executable(urdf-from-file apps/sample/trimesh_inertia/urdf_from_file.cpp wrap/ply/plylib.cpp)

install(TARGETS urdf-from-file DESTINATION ${CMAKE_INSTALL_PREFIX})

target_link_libraries(urdf-from-file Qt4::QtXml GLU)
