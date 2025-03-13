cmake_minimum_required(VERSION 3.8)

set(PROJECT_NAME centroidal_control_collection)
set(PROJECT_GENERATED_HEADERS_SKIP_DEPRECATED ON)
set(PROJECT_GENERATED_HEADERS_SKIP_CONFIG ON)
set(PROJECT_GENERATED_HEADERS_SKIP_WARNING ON)
set(PROJECT_URL https://github.com/isri-aist/CentroidalControlCollection)
set(PROJECT_DESCRIPTION "")
set(CMAKE_CXX_STANDARD 17)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CXX_DISABLE_WERROR ON)
set(CMAKE_COLOR_DIAGNOSTICS ON)

option(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)
option(USE_ROS2 "Use ROS2" OFF)
option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)

project(centroidal_control_collection LANGUAGES CXX)

if(USE_ROS2)
  # mc_rtc
  find_package(mc_rtc REQUIRED)

  find_package(rclcpp REQUIRED)
  find_package(qp_solver_collection REQUIRED)
  find_package(force_control_collection REQUIRED)
  find_package(nmpc_ddp REQUIRED)
else()
  set(BUILD_TESTING OFF)
  find_package(qp_solver_collection REQUIRED)
  find_package(force_control_collection REQUIRED)
  find_package(nmpc_ddp REQUIRED)
endif()

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(INSTALL_DOCUMENTATION)
  add_subdirectory(doc)
endif()

if(USE_ROS2)
  ament_export_include_directories(include)
  ament_export_libraries(CCC)
  ament_export_dependencies(
    rclcpp
    qp_solver_collection
    force_control_collection
    nmpc_ddp
  )
  ament_export_targets(CCC)
  ament_package()
endif() 