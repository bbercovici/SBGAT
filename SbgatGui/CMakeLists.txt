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
# @brief  CMake listing enabling compilation and installation of SBgatGui


################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################

cmake_minimum_required(VERSION 3.0.0)

# Building procedure
set(EXE_NAME "SbgatGui" CACHE STRING "Name of executable to be created.")

# Will use the most recent version of 
# Qt
SET(CMAKE_FIND_PACKAGE_SORT_ORDER NATURAL)
SET(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)

include(cmake/FindQt.cmake)

if(APPLE)
	set (CMAKE_PREFIX_PATH ${QT_CMAKE_CONFIG_MAC})
endif()

project(${EXE_NAME})

# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

add_definitions(-Wall -O2 )
include_directories(include)

find_package(Qt5 5.11.0 COMPONENTS Widgets REQUIRED)
include_directories(${Qt5Widgets_INCLUDE_DIRS})
add_definitions(${Qt5Widgets_DEFINITIONS})


# Find VTK Package
find_package(VTK REQUIRED )
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
include_directories(${RBK_INCLUDE_DIR})

# Find SBGAT and the JSON parsing library
find_package(SbgatCore REQUIRED)
include_directories(${SBGATCORE_INCLUDE_HEADER})
include_directories(${NLOHMANN_JSON_INCLUDE_HEADER})

# Find OrbitConversions 
find_package(OrbitConversions REQUIRED)
include_directories(${OC_INCLUDE_DIR})

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Add source files in root directory
add_executable(${EXE_NAME}
	source/main.cpp
	source/RenderingPropertiesWindow.cpp
	source/Mainwindow.cpp
	source/ModelDataWrapper.cpp
	source/SettingsWindow.cpp 
	source/YORPWindow.cpp
	source/SHARMWindow.cpp
	source/RadarWindow.cpp
	source/LCWindow.cpp
	source/RadarVisualizer.cpp
	source/LCVisualizer.cpp
	source/ObsWindow.cpp
	source/ShapePropertiesWidget.cpp
	source/SurfacePGMWindow.cpp
	source/SelectMapperWindow.cpp
	source/PickInteractorStyle.cpp
	source/VertexEditionWindow.cpp
	include/RenderingPropertiesWindow.hpp
	include/PickInteractorStyle.hpp
	include/VertexEditionWindow.hpp
	include/Mainwindow.hpp
	include/ModelDataWrapper.hpp
	include/SettingsWindow.hpp
	include/YORPWindow.hpp
	include/SHARMWindow.hpp
	include/RadarWindow.hpp
	include/LCWindow.hpp
	include/RadarVisualizer.hpp
	include/LCVisualizer.hpp
	include/ObsWindow.hpp
	include/ShapePropertiesWidget.hpp
	include/SurfacePGMWindow.hpp
	include/SelectMapperWindow.hpp

	)

set(library_dependencies
	${Qt5Widgets_LIBRARIES}
	${VTK_LIBRARIES}
	${ARMADILLO_LIBRARIES}
	${Boost_LIBRARIES}
	${SBGATCORE_LIBRARY}
	${RBK_LIBRARY}
	${OC_LIBRARY}
	)

target_link_libraries(${EXE_NAME} ${library_dependencies})


