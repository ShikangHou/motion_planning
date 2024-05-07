cmake_minimum_required(VERSION 3.0.2)
project(trajectory_generator)

## Compile as C++11, supported in ROS Kinetic and newer
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -O3 -march=native -g)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  visualization_msgs
)

find_package(Eigen3 REQUIRED)
set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

message(${EIGEN3_INCLUDE_DIR})

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES trajectory_generator
#  CATKIN_DEPENDS roscpp rospy std_msgs
#  DEPENDS system_lib
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS} 
)

add_executable(trajectory_generator_node src/trajectory_generator_node.cpp src/minimum_control.cpp)

target_link_libraries(trajectory_generator_node
                  ${catkin_LIBRARIES})