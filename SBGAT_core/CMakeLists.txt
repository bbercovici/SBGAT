# CMakeLists.txt for SBGAT_core
# Benjamin Bercovici, 03/10/2017
# ORCCA
# University of Colorado 


get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")
project(${EXE_NAME})

# Specify the version used


# cmake_minimum_required(VERSION 3.7.2)
message(CMAKE_VERSION)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)
add_definitions(-Wall)
include_directories(include)

if (UNIX AND NOT APPLE)
	set(CMAKE_CXX_FLAGS "-std=c++11 ")
elseif(APPLE)
	set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ")
	set(CMAKE_C_COMPILER clang-omp CACHE STRING "C compiler" FORCE)
	set(CMAKE_CXX_COMPILER clang-omp++ CACHE STRING "C++ compiler" FORCE)
endif()


# Find Boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS filesystem REQUIRED) 
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS}) 
else()
	message("Retrying find_boost")
	set(BOOST_ROOT "~/Libraries/boost/")
	set(BOOST_INCLUDEDIR "~/Libraries/boost/boost/")
	set(BOOST_LIBRARYDIR "~/Libraries/boost/libs/")


	find_package(Boost COMPONENTS filesystem REQUIRED) 
	include_directories(${Boost_INCLUDE_DIRS}) 
endif()


# Find Armadillo 
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Assimp
include_directories(~/include/assimp)
if (APPLE)
	set(ASSIMP_LIBRARY  
		"/usr/local/lib/libassimp.3.3.1.dylib"
		"/usr/local/lib/libassimp.3.dylib"
		"/usr/local/lib/libassimp.dylib")
elseif(UNIX AND NOT APPLE)
	set(ASSIMP_LIBRARY  
			"/usr/local/lib/libassimp.so.3.3.1"
			"/usr/local/lib/libassimp.so.3"
			"/usr/local/lib/libassimp.so")
endif()

# Fortran compiler, needed on Linux/Ubuntu by Armadillo
if(UNIX AND NOT APPLE)
	find_library(GFORTRAN_LIBRARY gfortran
	    PATHS /usr/lib/gcc/x86_64-linux-gnu/5/)
	list(APPEND ARMADILLO_LIBRARIES "${GFORTRAN_LIBRARY}")
endif()

# Add source files in root directory
add_executable(${EXE_NAME}
source/main.cpp
source/ShapeModel.cpp
source/DynamicAnalyses.cpp
include/ShapeModel.hpp
include/DynamicAnalyses.hpp
include/Densities.hpp
)

set(library_dependencies
${ARMADILLO_LIBRARIES}
${ASSIMP_LIBRARY}
${Boost_LIBRARIES})

target_link_libraries(${EXE_NAME} ${library_dependencies})

