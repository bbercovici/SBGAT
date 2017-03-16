# CMakeLists.txt for SBGAT_core
# Benjamin Bercovici, 03/10/2017
# ORCCA
# University of Colorado 


get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")


set(GOLEVKA_NAME Linux-2.6.32-642.11.1.el6.x86_64)
message("${CMAKE_SYSTEM}")

if ("${CMAKE_SYSTEM}" STREQUAL "${GOLEVKA_NAME}")
	set(CMAKE_C_COMPILER /usr/local/gcc/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/local/gcc/bin/g++ CACHE STRING "C++ Compiler" FORCE)
else()
	message("${CMAKE_SYSTEM}")
	message("${GOLEVKA_NAME}")
	message("This is not Golevka")
endif()


project(${EXE_NAME})

message(${CMAKE_SYSTEM})


# Specify the version used

if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
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



# Check what version of the compiler is installed
execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                OUTPUT_VARIABLE GCC_VERSION)
message(${GCC_VERSION})
if (GCC_VERSION VERSION_LESS 4.7)
        set(CMAKE_CXX_FLAGS "-std=c++0x")
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
	set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /homes/bebe0705/Libraries/boost)
	set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /homes/bebe0705/Libraries/boost/libs/)
	message( ${CMAKE_LIBRARY_PATH})
	message( ${CMAKE_INCLUDE_PATH})
	
	find_package(Boost COMPONENTS filesystem REQUIRED PATHS ~/Libraries/boost/) 



# 	set(Boost_INCLUDE_DIR /usr/local/src/boost_1_46_1)
# set(Boost_LIBRARY_DIR /usr/local/src/boost_1_46_1/stage/lib)
# find_package(Boost COMPONENTS system filesystem REQUIRED)
# include_directories(${Boost_INCLUDE_DIR})
# link_directories(${Boost_LIBRARY_DIR})



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
# find_package(assimp)
# include_directories(${ASSIMP_INCLUDE_DIRS})


# include_directories(~/include/assimp)
# if (APPLE)
# 	set(ASSIMP_LIBRARY  
# 		"/usr/local/lib/libassimp.3.3.1.dylib"
# 		"/usr/local/lib/libassimp.3.dylib"
# 		"/usr/local/lib/libassimp.dylib")
# elseif(UNIX AND NOT APPLE)
# 	set(ASSIMP_LIBRARY  
# 			"/usr/local/lib/libassimp.so.3.3.1"
# 			"/usr/local/lib/libassimp.so.3"
# 			"/usr/local/lib/libassimp.so")
# endif()

# Fortran compiler, needed on Linux/Ubuntu by Armadillo
if(UNIX AND NOT APPLE)
	find_library(GFORTRAN_LIBRARY gfortran
	    PATHS /usr/lib/gcc/x86_64-linux-gnu/5/ /usr/lib/gcc/x86_64-redhat-linux/4.4.7/32/)
	list(APPEND ARMADILLO_LIBRARIES "${GFORTRAN_LIBRARY}")
endif()

# Add source files in root directory
add_executable(${EXE_NAME}
source/main.cpp
source/ShapeModel.cpp
source/DynamicAnalyses.cpp
source/ShapeModelImporter.cpp
include/ShapeModel.hpp
include/DynamicAnalyses.hpp
include/Densities.hpp
include/ShapeModelImporter.hpp
)

set(library_dependencies
${ARMADILLO_LIBRARIES}
${Boost_LIBRARIES})

target_link_libraries(${EXE_NAME} ${library_dependencies})

