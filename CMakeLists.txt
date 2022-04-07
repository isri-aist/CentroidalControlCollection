cmake_minimum_required(VERSION 3.1)
project(centroidal_control_collection)

add_compile_options(-std=c++17)

find_package(catkin REQUIRED COMPONENTS
  message_generation
  roscpp
  qp_solver_collection
  nmpc_ddp
  )

catkin_package(
  CATKIN_DEPENDS
  roscpp
  qp_solver_collection
  nmpc_ddp
  DEPENDS EIGEN3
  INCLUDE_DIRS include
  LIBRARIES CCC
  )

add_subdirectory(src)

add_subdirectory(samples)

if(CATKIN_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

OPTION(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)
if(INSTALL_DOCUMENTATION)
  add_subdirectory(doc)
endif()
