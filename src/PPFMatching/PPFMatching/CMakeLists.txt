cmake_minimum_required (VERSION 2.8)
project(PPFMatching)

find_package(OpenCV REQUIRED)
find_package(Eigen REQUIRED)
#find_package(flann REQUIRED)

#include_directories(${OpenCV_INCLUDE_DIR})
include_directories(${EIGEN_INCLUDE_DIRS} 
                    includes)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_library(matrix3x3
  dsyevh3.c
  dsyevq3.c
  dsyevc3.c
  dsytrd3.c
)

add_executable(main 
  main.cpp
  helpers.cpp
  t_octree.cpp
  tommyhashtbl.c
  fasthash.c
  t_icp.cpp
  THashInt.c
  c_utils.cpp
  gdiam.cpp
)
target_link_libraries(main ${OpenCV_LIBS}
  flann matrix3x3)
