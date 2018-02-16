# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   January, 2018
# @brief  CMake listing enabling compilation of PolyhedronGravityModel example


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
set (OMP_FRIENDLY_GCC_PATH /usr/local/Cellar/gcc/7.2.0/bin/)


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
	set(RBK_LOC "/home/bebe0705/libs/local/lib/cmake/RigidBodyKinematics")
	set(SBGAT_LOC "/home/bebe0705/libs/local/lib/cmake/SbgatCore")

else()
	set(IS_FORTUNA OFF)
	message("-- This is not Fortuna")

endif()


# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")


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



project(${EXE_NAME})

# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)


add_definitions(-Wall -O2 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

include_directories(include)

# Find Boost
find_package(Boost COMPONENTS filesystem) 
include_directories(${Boost_INCLUDE_DIRS}) 

# Find armadillo package
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find RBK 
find_package(RigidBodyKinematics REQUIRED PATHS ${RBK_LOC})
include_directories(${RBK_INCLUDE_DIR})

# Find SBGAT 
find_package(SbgatCore REQUIRED PATHS ${SBGAT_LOC})
include_directories(${SBGATCORE_INCLUDE_HEADER})


# Find OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


# Add source files in root directory
add_executable(${EXE_NAME}
source/main.cpp
)

set(library_dependencies
${ARMADILLO_LIBRARIES}
${Boost_LIBRARIES}
${SBGATCORE_LIBRARY}
${RBK_LIBRARY}
)

target_link_libraries(${EXE_NAME} ${library_dependencies})

