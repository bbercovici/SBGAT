get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")
project(${EXE_NAME})

# Specify the version used
cmake_minimum_required(VERSION 3.0.2)


set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)


add_definitions(-Wall)
include_directories(include)



# Find Boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS Filesystem) 
include_directories(${Boost_INCLUDE_DIRS}) 


# Find armadillo package
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Assimp
include_directories(/Users/bbercovici/assimp_3_3_1/include)
set(assimp_dependencies  
	"/usr/local/lib/libassimp.3.3.1.dylib"
	"/usr/local/lib/libassimp.3.dylib"
	"/usr/local/lib/libassimp.dylib")




# Add source files in root directory
add_executable(${EXE_NAME}
source/main.cpp
source/ShapeModel.cpp
source/DynamicAnalyses.cpp
include/ShapeModel.hpp
include/DynamicAnalyses.hpp
)

set(library_dependencies
${ARMADILLO_LIBRARIES}
${assimp_dependencies}
${Boost_LIBRARIES})


target_link_libraries(${EXE_NAME} ${library_dependencies})

