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


# OMP-friendly GCC Major version (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_MAJOR_VERSION 7)

# Path to OMP-Friendly GCC Compiler (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_PATH /usr/local/Cellar/gcc/7.3.0_1/bin/)

# Path to VTKConfig.cmake (not mandatory if VTK is installed at a standard location)
set (VTK_PATH /usr/local/VTK-8.1.0/lib/cmake/vtk-8.1)


################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################


if (EXISTS /home/bebe0705/.am_fortuna)
	set(IS_FORTUNA ON)
	message("-- This is Fortuna")

else()
	set(IS_FORTUNA OFF)
	message("-- This is not Fortuna")

endif()



# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME ${dirName} CACHE STRING "Name of executable to be created.")


if (${IS_FORTUNA})
	set(CMAKE_C_COMPILER "/usr/local/bin/gcc" CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER "/usr/local/bin/g++" CACHE STRING "C++ Compiler" FORCE)
else()
	if(APPLE)
		if(EXISTS ${OMP_FRIENDLY_GCC_PATH}) # Else running on a MAC. Will attempt to switch compiler to get the OMP-friendly GCC 6.3.0_1 from Homebrew
			message("Switching to OMP-friendly GCC ")
	    	set(CMAKE_C_COMPILER ${OMP_FRIENDLY_GCC_PATH}gcc-${OMP_FRIENDLY_GCC_MAJOR_VERSION} CACHE STRING "C Compiler" FORCE)
			set(CMAKE_CXX_COMPILER ${OMP_FRIENDLY_GCC_PATH}g++-${OMP_FRIENDLY_GCC_MAJOR_VERSION} CACHE STRING "C++ Compiler" FORCE)

		else()# no OMP-compliant compiler was found on this mac. 
			message(FATAL_ERROR "No OMP-compliant compiler was found on this Mac.")
		endif()
	else() # Running on Linux. Will switch back to compiler in /usr/local/bin
				message("Switching to /usr/local/gcc ")
		    	set(CMAKE_C_COMPILER "/usr/local/bin/gcc" CACHE STRING "C Compiler" FORCE)
				set(CMAKE_CXX_COMPILER "/usr/local/bin/g++" CACHE STRING "C++ Compiler" FORCE)
	endif()
endif()

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# Include directories
include_directories(include)

# Find Boost
find_package(Boost COMPONENTS filesystem system REQUIRED) 
include_directories(${Boost_INCLUDE_DIRS}) 

# Find Armadillo 
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Find VTK Package
find_package(VTK REQUIRED PATHS ${VTK_PATH})
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
if(UNIX AND NOT APPLE AND ${CMAKE_MINOR_VERSION} GREATER 0 AND NOT ${IS_FORTUNA})
	find_library(GFORTRAN_LIBRARY gfortran
	    PATHS /usr/lib/gcc/x86_64-linux-gnu/5/ /usr/lib/gcc/x86_64-redhat-linux/4.4.7/32/)
	list(APPEND ARMADILLO_LIBRARIES "${GFORTRAN_LIBRARY}")
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
if (NOT ${IS_FORTUNA})
	install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib/)
	install (FILES ${PROJECT_SOURCE_DIR}/SbgatCoreConfig.cmake DESTINATION /usr/local/lib/cmake/SbgatCore/)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION /usr/local/include/SbgatCore FILES_MATCHING PATTERN "*.hpp")
else()
	install (TARGETS ${LIB_NAME} DESTINATION /home/bebe0705/libs/local/lib)
	install (FILES ${PROJECT_SOURCE_DIR}/SbgatCoreConfig.cmake DESTINATION /home/bebe0705/libs/local/lib/cmake/SbgatCore/)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION /home/bebe0705/libs/local/include/SbgatCore FILES_MATCHING PATTERN "*.hpp")
endif()







