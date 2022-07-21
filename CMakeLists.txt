cmake_minimum_required(VERSION 2.9 FATAL_ERROR)

project(mesh LANGUAGES CXX )

MESSAGE("CMAKE_CXX_COMPILER_VERSION: " ${CMAKE_CXX_COMPILER_VERSION})


add_library(extras_lib STATIC
	TissueBacteria.cpp
	Bacteria.cpp
	Nodes.cpp
	SignalCommon.cpp
	ConfigParser.cpp
	commonData.cpp
	GeoVector.cpp
	ranum2.cpp
	Diffusion2D.cpp
	Grid.cpp
	TissueGrid.cpp
	driver.cpp
	Fungi.cpp
	growthFunctions.cpp
	hyphaeSegment.cpp
	constants.cpp)

	
target_compile_features(extras_lib PUBLIC cxx_std_11)

set(MAIN_SOURCES main.cpp)
set(TARGET program)
add_executable(${TARGET} ${MAIN_SOURCES})

MESSAGE("TARGET: " ${TARGET})

target_link_libraries(${TARGET} extras_lib )

