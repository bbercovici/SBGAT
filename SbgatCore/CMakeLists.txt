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

# Path to OMP-Friendly GCC Compiler (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_PATH /usr/local/Cellar/gcc/6.3.0_1/bin/)




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

if (${CMAKE_MINOR_VERSION} EQUAL 0 AND NOT APPLE) # if true, running on a CCAR computer (old CMAKE and RedHat)
	set(CMAKE_C_COMPILER /usr/local/gcc/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/local/gcc/bin/g++ CACHE STRING "C++ Compiler" FORCE)
elseif(NOT APPLE)
	# Not sure about this. Should work on a modern Ubuntu system
	set(CMAKE_C_COMPILER /usr/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/bin/g++ CACHE STRING "C++ Compiler" FORCE)
    set (CMAKE_CXX_STANDARD 17)
	message("This is a Linux machine. CMAKE_CXX_STANDARD 17 used")


endif()

# Check what version of the compiler is installed
execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                OUTPUT_VARIABLE GCC_VERSION)
message("Default compiler is " ${CMAKE_C_COMPILER})
message("Default ${CMAKE_C_COMPILER} version is " ${GCC_VERSION})

# Set up the C++11 and OMP compile flags 
if (GCC_VERSION VERSION_LESS 4.7 OR APPLE) # if true, either running on a CCAR computer or a Mac (Clang has messed-up version IDs)
		if (NOT APPLE) # If NOT APPLE, then running on a CCAR computer. The proper compiler should have already been set
        	set(CMAKE_CXX_FLAGS "-std=c++0x")
        elseif(EXISTS ${OMP_FRIENDLY_GCC_PATH}) # Else running on a MAC. Will attempt to switch compiler to get the OMP-friendly GCC 6.3.0_1 from Homebrew
        		message("Running on Mac ")
        		message("Switching to OMP-friendly GCC ")
	        	set(CMAKE_C_COMPILER ${OMP_FRIENDLY_GCC_PATH}gcc-6 CACHE STRING "C Compiler" FORCE)
				set(CMAKE_CXX_COMPILER ${OMP_FRIENDLY_GCC_PATH}g++-6 CACHE STRING "C++ Compiler" FORCE)
			else()
				message(FATAL_ERROR "No OMP-compliant compiler was found on this mac. ")
		endif()
			set(CMAKE_CXX_FLAGS "-std=c++11")

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


# Find RBK 
# set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} ${RBK_LIBRARY_CMAKE_INSTALL}")
find_package(RigidBodyKinematics REQUIRED)
include_directories(${RBK_INCLUDE_DIRS})

# Fortran compiler, required by Armadillo on Linux/Ubuntu
if(UNIX AND NOT APPLE AND ${CMAKE_MINOR_VERSION} GREATER 0)
	find_library(GFORTRAN_LIBRARY gfortran
	    PATHS /usr/lib/gcc/x86_64-linux-gnu/5/ /usr/lib/gcc/x86_64-redhat-linux/4.4.7/32/)
	list(APPEND ARMADILLO_LIBRARIES "${GFORTRAN_LIBRARY}")
endif()

# Add source files in root directory
add_library(${LIB_NAME} 
SHARED
source/ShapeModel.cpp
source/DynamicAnalyses.cpp
source/ShapeModelImporter.cpp
source/Facet.cpp
source/Edge.cpp
source/Vertex.cpp
source/RefFrame.cpp
source/FrameGraph.cpp
source/SyntheticObservations.cpp
source/FacetResults.cpp
)

# Linking
set(library_dependencies
${ARMADILLO_LIBRARIES}
${Boost_LIBRARIES}
${RBK_LIBRARY})
target_link_libraries(${LIB_NAME} ${library_dependencies})

# Installing
install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib/)
install (FILES ${PROJECT_SOURCE_DIR}/SbgatCoreConfig.cmake DESTINATION /usr/local/lib/cmake/SbgatCore/)
install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION /usr/local/include/SbgatCore FILES_MATCHING PATTERN "*.hpp")








