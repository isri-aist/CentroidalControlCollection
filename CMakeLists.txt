cmake_minimum_required(VERSION 3.8)

set(PROJECT_NAME centroidal_control_collection)
set(PROJECT_GENERATED_HEADERS_SKIP_DEPRECATED ON)
set(PROJECT_GENERATED_HEADERS_SKIP_CONFIG ON)
set(PROJECT_GENERATED_HEADERS_SKIP_WARNING ON)
set(PROJECT_URL https://github.com/isri-aist/centroidal_control_collection)
set(PROJECT_DESCRIPTION "")
set(CMAKE_CXX_STANDARD 17)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CXX_DISABLE_WERROR ON)
option(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)

include(cmake/base.cmake)
project(centroidal_control_collection LANGUAGES CXX)

if(DEFINED CATKIN_DEVEL_PREFIX)
  set(DOXYGEN_HTML_OUTPUT html)
  find_package(catkin REQUIRED COMPONENTS
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
else()
  set(DOXYGEN_HTML_OUTPUT doxygen-html)
  set(CATKIN_ENABLE_TESTING OFF)
  option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)
  add_project_dependency(qp_solver_collection REQUIRED)
  add_project_dependency(nmpc_ddp REQUIRED)
endif()

add_subdirectory(src)

add_subdirectory(samples)

if(BUILD_TESTING OR CATKIN_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

OPTION(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)
if(INSTALL_DOCUMENTATION)
  add_subdirectory(doc)
endif()
