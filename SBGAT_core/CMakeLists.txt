# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   July, 2017
# @brief  CMake listing enabling compilation and installation of the SBGAT_CORE library


# Specify install location by setting this variable. Currently installing under
# ./SBGAT_core/lib/
set (SBGAT_CORE_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib )

# Specify location of RBK Library headers
set(RBK_LIBRARY_CMAKE_INSTALL /usr/local/lib/cmake/RigidBodyKinematics/)


# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME ${dirName} CACHE STRING "Name of executable to be created.")

if (${CMAKE_MINOR_VERSION} EQUAL 0 AND NOT APPLE) # if true, running on a CCAR computer (old CMAKE and RedHat)
	set(CMAKE_C_COMPILER /usr/local/gcc/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/local/gcc/bin/g++ CACHE STRING "C++ Compiler" FORCE)
else()
	# Not sure about this. Should work on a modern Ubuntu system
	set(CMAKE_C_COMPILER /usr/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/bin/g++ CACHE STRING "C++ Compiler" FORCE)
	message("This is not a CCAR computer")
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
        else() # Else running on a MAC. Will attempt to switch compiler to get the OMP-friendly GCC 6.3.0_1 from Homebrew
        	if(EXISTS /usr/local/Cellar/gcc/6.3.0_1/bin/) 
        		message("Switching to GCC 6.3.0_1")
	        	set(CMAKE_C_COMPILER /usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6 CACHE STRING "C Compiler" FORCE)
				set(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/6.3.0_1/bin/g++-6 CACHE STRING "C++ Compiler" FORCE)
			elseif(EXISTS /usr/local/bin/clang-omp) # If the previous check failed, will attempt to get clang-omp/clang-omp++
			    message("Switching to Clamg-OMP")
				set(CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C Compiler" FORCE)
				set(CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ Compiler" FORCE)
			else()# no OMP-compliant compiler was found on this mac. The program will still run but without multithreading
				message("No OMP-compliant compiler was found on this mac. The program will still run but without multithreading")
			endif()
			set(CMAKE_CXX_FLAGS "-std=c++11")

        endif()
else() # Running on a proper Linux machine
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
include_directories(external/include)

# Find Boost
set(BOOST_ROOT /usr/local/Cellar/boost/1.63.0/)
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
set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} ${RBK_LIBRARY_CMAKE_INSTALL}")
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
)

# Linking
set(library_dependencies
${ARMADILLO_LIBRARIES}
${Boost_LIBRARIES}
${RBK_LIBRARY})
target_link_libraries(${LIB_NAME} ${library_dependencies})

# Installing
install (TARGETS ${LIB_NAME} DESTINATION ${SBGAT_CORE_INSTALL_DIR})
install (FILES ${RBK_LIBRARY_CMAKE_INSTALL}/RigidBodyKinematicsConfig.cmake DESTINATION ${SBGAT_CORE_INSTALL_DIR})
install (FILES ${PROJECT_SOURCE_DIR}/SBGATCoreConfig.cmake DESTINATION /usr/local/lib/cmake/SBGAT_CORE/)








