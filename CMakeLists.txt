cmake_minimum_required(VERSION 3.1)
project(centroidal_control_collection)

add_compile_options(-std=c++17)

find_package(catkin REQUIRED COMPONENTS
  message_generation
  roscpp
  rosbag
  rviz
  urdf
  std_msgs
  std_srvs
  eigen_conversions
  geometry_msgs
  visualization_msgs
  sensor_msgs
  qp_solver_collection
  )

catkin_package(
  CATKIN_DEPENDS
  message_runtime
  roscpp
  rosbag
  rviz
  urdf
  std_msgs
  std_srvs
  eigen_conversions
  geometry_msgs
  visualization_msgs
  sensor_msgs
  qp_solver_collection
  DEPENDS EIGEN3
  INCLUDE_DIRS include
  LIBRARIES CCC
  )

add_subdirectory(src)

add_subdirectory(samples)

if(CATKIN_ENABLE_TESTING)
  add_subdirectory(tests)
endif()
