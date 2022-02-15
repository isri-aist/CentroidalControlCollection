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
  optmotiongen
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
  optmotiongen
  DEPENDS EIGEN3
  LIBRARIES CCC
  )

include_directories(include ${catkin_INCLUDE_DIRS})
link_directories(${catkin_LIBRARY_DIRS})

add_subdirectory(src)

add_subdirectory(tests)
