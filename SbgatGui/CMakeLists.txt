



get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")
project(${EXE_NAME})

# Specify the version used
cmake_minimum_required(VERSION 3.0.2)

# Set project name - based on current directory
# Should give you the option to change in the CMake gui

# Specify the available build types for multi-configuration generators
# set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE INTERNAL "Configuration types")

# Specify the directory that simulation code is located in
# Should give you the option to change in the CMake gui
# set(SIMULATION_DIR ${CMAKE_SOURCE_DIR}/../Simulation/ CACHE PATH "Directory containing simulation source code")
# Put all the runtime targets in one folder
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin/Debug")
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin/Release")

set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

if (UNIX AND NOT APPLE)
	set (CMAKE_PREFIX_PATH "~/Qt5.7.1/5.7/gcc_64/lib/cmake/" "~/VTK/VTK_release_build/")
	set(CMAKE_CXX_FLAGS "-std=c++11 ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

elseif(APPLE)
	set (CMAKE_PREFIX_PATH "~/Qt/5.9.1/clang_64/lib/cmake" "~/VTK/VTK_release_build/")
	set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

endif()


add_definitions(-Wall -O2)
include_directories(include)

# Widgets finds its own dependencies (QtGui and QtCore)
find_package(Qt5Widgets REQUIRED)

# Find VTK Package
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

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

# Cocoa
if(APPLE)
	find_library( LIB_COCOA cocoa )
endif()

# The Qt5Widgets_INCLUDES also includes the include directories for dependencies QtCore and QtGui
include_directories(${Qt5Widgets_INCLUDES})

# We need add -DQT_WIDGETS_LIB when using QtWidgets in Qt 5.
add_definitions(${Qt5Widgets_DEFINITIONS})

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)# As moc files are generated in the binary dir, tell CMake to always look for includes there

# Add source files in root directory
add_executable(${EXE_NAME}
source/main.cpp
source/Mainwindow.cpp
source/osxHelper.mm
include/Mainwindow.hpp
include/osxHelper.h   
)

set(library_dependencies
${Qt5Widgets_LIBRARIES}
${VTK_LIBRARIES}
${ARMADILLO_LIBRARIES}
${Boost_LIBRARIES})


target_link_libraries(${EXE_NAME} ${library_dependencies})

