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
# @brief  CMake listing enabling compilation and installation of the SbgatCore tests

################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################

set (RBK_PATH "/home/bebe0705/libs/local/lib/cmake/RigidBodyKinematics")
set (SBGAT_PATH "/home/bebe0705/libs/local/lib/cmake/SbgatCore")
set (VTK_PATH "/usr/local/VTK-8.1.0/lib/cmake/vtk-8.1")



cmake_minimum_required(VERSION 3.05.0)

# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")


project(${EXE_NAME})

# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

add_definitions(-Wall -O2 )
include_directories(include)

# Find Boost
find_package(Boost COMPONENTS filesystem) 
include_directories(${Boost_INCLUDE_DIRS}) 

# Find armadillo package
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find RBK 
find_package(RigidBodyKinematics REQUIRED PATHS ${RBK_PATH})
include_directories(${RBK_INCLUDE_DIR})

# Find SBGAT 
find_package(SbgatCore REQUIRED PATHS ${SBGAT_PATH})
include_directories(${SBGATCORE_INCLUDE_HEADER})
include_directories(${NLOHMANN_JSON_INCLUDE_HEADER})
# Find OpenMP
if (NOT NO_OMP)
	find_package(OpenMP)

	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()

endif()



# Find VTK Package
find_package(VTK REQUIRED PATHS ${VTK_PATH})
include(${VTK_USE_FILE})

# Removing spurious include sometimes brought in by one of VTK's dependencies
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
list(REMOVE_ITEM dirs "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${dirs})

# Add source files in root directory
add_executable(${EXE_NAME} main.cpp )


set(library_dependencies
	${ARMADILLO_LIBRARIES}
	${Boost_LIBRARIES}
	${SBGATCORE_LIBRARY}
	${RBK_LIBRARY}
	${VTK_LIBRARIES}
	)



target_link_libraries(${EXE_NAME} ${library_dependencies})


