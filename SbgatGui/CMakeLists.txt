# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   July, 2017
# @brief  CMake listing enabling compilation and installation of SBgatGui


################################################################################
#
# 								User-defined paths
#						Should be checked for consistency
#						Before running 'cmake ..' in build dir
#
################################################################################

# Path to OMP-Friendly GCC Compiler (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_PATH /usr/local/Cellar/gcc/6.3.0_1/bin/)

# Path to VTK build directory
set (VTK_BUILD_DIR ~/VTK/VTK_release_build/)

# Path to Qt5 cmake config file (author's mac configuration)
set (QT_CMAKE_CONFIG_MAC "~/Qt/5.9.1/clang_64/lib/cmake")

# Path to Qt5 cmake config file (author's ubuntu configuration)
set (QT_CMAKE_CONFIG_LINUX "~/Qt/5.9.1/gcc_64/lib/cmake/")



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
set(EXE_NAME ${dirName} CACHE STRING "Name of executable to be created.")

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
        		message("Switching to OMP-friendly GCC ")
	        	set(CMAKE_C_COMPILER ${OMP_FRIENDLY_GCC_PATH}gcc-6 CACHE STRING "C Compiler" FORCE)
				set(CMAKE_CXX_COMPILER ${OMP_FRIENDLY_GCC_PATH}g++-6 CACHE STRING "C++ Compiler" FORCE)
			else()# no OMP-compliant compiler was found on this mac. 
				message(FATAL_ERROR "No OMP-compliant compiler was found on this mac. The program will still run but without multithreading")
		endif()
		set(CMAKE_CXX_FLAGS "-std=c++17")

endif()


project(${EXE_NAME})

# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)


if (UNIX AND NOT APPLE)
	set (CMAKE_PREFIX_PATH ${QT_CMAKE_CONFIG_LINUX} ${VTK_BUILD_DIR})
	set(CMAKE_CXX_FLAGS "-std=c++11 ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

elseif(APPLE)
	set (CMAKE_PREFIX_PATH ${QT_CMAKE_CONFIG_MAC} ${VTK_BUILD_DIR})
	set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
endif()


add_definitions(-Wall -O2 )
set(CMAKE_CXX_FLAGS "-std=c++11")
include_directories(include)

# Widgets finds its own dependencies (QtGui and QtCore)
find_package(Qt5Widgets REQUIRED)

# The Qt5Widgets_INCLUDES also includes the include directories for dependencies QtCore and QtGui
include_directories(${Qt5Widgets_INCLUDES})

# We need add -DQT_WIDGETS_LIB when using QtWidgets in Qt 5.
add_definitions(${Qt5Widgets_DEFINITIONS})

# Find VTK Package
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Find Boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS filesystem) 
include_directories(${Boost_INCLUDE_DIRS}) 


# Find armadillo package
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find RBK 
find_package(RigidBodyKinematics REQUIRED)
include_directories(${RBK_INCLUDE_DIRS})

# Find SBGAT 
find_package(SbgatCore REQUIRED)
include_directories(${SBGATCORE_INCLUDE_HEADER})

# Find OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)# As moc files are generated in the binary dir, tell CMake to always look for includes there

# Add source files in root directory
add_executable(${EXE_NAME}
source/main.cpp
source/Mainwindow.cpp
source/SettingsWindow.cpp
source/ModelDataWrapper.cpp
source/Worker.cpp
include/Mainwindow.hpp
include/SettingsWindow.hpp
include/ModelDataWrapper.hpp
include/Worker.hpp
)

set(library_dependencies
${Qt5Widgets_LIBRARIES}
${VTK_LIBRARIES}
${ARMADILLO_LIBRARIES}
${Boost_LIBRARIES}
${SBGATCORE_LIBRARY}
${RBK_LIBRARY}
)

target_link_libraries(${EXE_NAME} ${library_dependencies})

