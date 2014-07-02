cmake_minimum_required (VERSION 2.8)
project(PPFMatching)

find_package(OpenCV REQUIRED)
find_package(Eigen REQUIRED)

#include_directories(${OpenCV_INCLUDE_DIR})
include_directories(${EIGEN_INCLUDE_DIRS} 
					includes)

add_executable(main 
	main.cpp
	helpers.cpp
	t_octree.cpp
	tommyhashtbl.c
	gdiam.cpp
	c_utils.cpp
	)
target_link_libraries(main ${OpenCV_LIBS})
