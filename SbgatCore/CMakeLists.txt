# MIT License

# Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   July, 2017
# @brief  CMake listing enabling compilation and installation of the SbgatCore library


################################################################################
#
# 								User-defined paths
#						Should be checked for consistency
#						Before running 'cmake ..' in build dir
#
################################################################################


################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################
	
cmake_minimum_required(VERSION 3.0.0)

if (${USE_GCC})
	include(cmake/FindOmpGcc.cmake)
else()
	set(CMAKE_C_COMPILER /usr/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/bin/g++ CACHE STRING "C++ Compiler" FORCE)
endif()

# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME "SbgatCore" CACHE STRING "Name of executable to be created.")

project(${LIB_NAME})

# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

# Compiler flags
add_definitions(-Wall -O2 )

# Enable C++17 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Include directories
include_directories(include/SbgatCore)
include_directories(include/nlohmann)

# Find Boost
find_package(Boost COMPONENTS filesystem system REQUIRED) 
include_directories(${Boost_INCLUDE_DIRS}) 

# Find Armadillo 
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Find VTK Package
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Find YORPLib
find_package(YORPLib REQUIRED)
include_directories(${YORPLIB_INCLUDE_HEADER})

# Find SHARMLib
find_package(SHARMLib REQUIRED)
include_directories(${SHARMLIB_INCLUDE_HEADER})

# Find RBK 
find_package(RigidBodyKinematics REQUIRED )
include_directories(${RBK_INCLUDE_DIR})

# Fortran compiler, required by Armadillo on Linux/Ubuntu
if(UNIX AND NOT APPLE AND ${CMAKE_MINOR_VERSION} GREATER 0 )
	find_library(GFORTRAN_LIBRARY gfortran
		PATHS /usr/lib/gcc/x86_64-linux-gnu/5/ /usr/lib/gcc/x86_64-redhat-linux/4.4.7/32/)
	if (GFORTRAN_LIBRARY_FOUND)
		list(APPEND ARMADILLO_LIBRARIES "${GFORTRAN_LIBRARY}")
	endif()
endif()

# Add source files in root directory
add_library(${LIB_NAME} 
	SHARED
	source/SBGATRefFrame.cpp
	source/SBGATFrameGraph.cpp
	source/SBGATPolyhedronGravityModel.cpp
	source/SBGATMassProperties.cpp
	source/SBGATSrpYorp.cpp
	source/SBGATSphericalHarmo.cpp
	source/SBGATObjWriter.cpp
	source/SBGATObsRadar.cpp
	source/SBGATObsLightcurve.cpp

	)

# Linking
set(library_dependencies
	${ARMADILLO_LIBRARIES}
	${Boost_LIBRARIES}
	${VTK_LIBRARIES}
	${RBK_LIBRARY}
	${YORPLIB_LIBRARY}
	${SHARMLIB_LIBRARY}
	)
target_link_libraries(${LIB_NAME} ${library_dependencies})

# Installing
if(NOT BREW)
	install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib/)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/${LIB_NAME} DESTINATION /usr/local/share)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/SbgatCore DESTINATION /usr/local/include FILES_MATCHING PATTERN "*.hpp")
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/nlohmann DESTINATION /usr/local/include FILES_MATCHING PATTERN "*.hpp")
endif()







