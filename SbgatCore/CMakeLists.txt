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

# Path to VTKConfig.cmake (not mandatory if VTK is installed at a standard location)
set (VTK_PATH /usr/local/VTK-8.1.0/lib/cmake/vtk-8.1)
set(RBK_PATH "/home/bebe0705/libs/local/lib/cmake/RigidBodyKinematics")
set (YORPLIB_PATH "/home/bebe0705/libs/local/lib/cmake/YORPLib")
set (SHARMLIB_PATH "/home/bebe0705/libs/local/lib/cmake/SHARMLib")
set (NLOHMANN_JSON_PATH "/home/bebe0705/libs/local/lib/cmake/nlohmann_json")



################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################

# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME ${dirName} CACHE STRING "Name of executable to be created.")

set(CMAKE_C_COMPILER "/usr/local/bin/gcc" CACHE STRING "C Compiler" FORCE)
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++" CACHE STRING "C++ Compiler" FORCE)

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
include_directories(include/SbgatCore)
include_directories(include/nlohmann)

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
find_package(YORPLib REQUIRED PATHS ${YORPLIB_PATH})
include_directories(${YORPLIB_INCLUDE_HEADER})


# Find OrbitConversions 
find_package(OrbitConversions REQUIRED)
include_directories(${OC_INCLUDE_DIR})

# Find SHARMLib
find_package(SHARMLib REQUIRED PATHS ${SHARMLIB_PATH})
include_directories(${SHARMLIB_INCLUDE_HEADER})

# Find nlohmann_json
find_package(nlohmann_json REQUIRED PATHS ${NLOHMANN_JSON_PATH})

# Find RBK 
find_package(RigidBodyKinematics REQUIRED PATHS ${RBK_PATH})
include_directories(${RBK_INCLUDE_DIR})

# Add source files in root directory
add_library(${LIB_NAME} 
	SHARED
	source/SBGATRefFrame.cpp
	source/SBGATFilter.cpp
	source/SBGATFrameGraph.cpp
	source/SBGATPolyhedronGravityModel.cpp
	source/SBGATPolyhedronGravityModelUQ.cpp
	source/SBGATMassProperties.cpp
	source/SBGATMassPropertiesUQ.cpp
	source/SBGATShapeUncertainty.cpp
	source/SBGATSrpYorp.cpp
	source/SBGATSphericalHarmo.cpp
	source/SBGATObjWriter.cpp
	source/SBGATObs.cpp
	source/SBGATObsRadar.cpp
	source/SBGATObsLightcurve.cpp
	source/SBGATTrajectory.cpp
	source/SBGATTransformShape.cpp)

# Linking
set(library_dependencies
	${ARMADILLO_LIBRARIES}
	${Boost_LIBRARIES}
	${VTK_LIBRARIES}
	${RBK_LIBRARY}
	${OC_LIBRARY}
	${YORPLIB_LIBRARY}
	${SHARMLIB_LIBRARY}
	)
target_link_libraries(${LIB_NAME} ${library_dependencies})

# Installing
install (TARGETS ${LIB_NAME} DESTINATION /home/bebe0705/libs/local/lib)
# install (FILES ${PROJECT_SOURCE_DIR}/SbgatCore/CMake/SbgatCoreConfig.cmake DESTINATION /home/bebe0705/libs/local/lib/cmake/SbgatCore/)
# install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION /home/bebe0705/libs/local/include/SbgatCore FILES_MATCHING PATTERN "*.hpp")
install (DIRECTORY ${PROJECT_SOURCE_DIR}/${LIB_NAME} DESTINATION /home/bebe0705/libs/local/lib/cmake)
install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/SbgatCore DESTINATION /home/bebe0705/libs/local/include FILES_MATCHING PATTERN "*.hpp")
install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/nlohmann DESTINATION /home/bebe0705/libs/local/include FILES_MATCHING PATTERN "*.hpp")




