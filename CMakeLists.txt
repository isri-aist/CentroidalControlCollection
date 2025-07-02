cmake_minimum_required(VERSION 3.10)

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
include(GNUInstallDirs) # For CMAKE_INSTALL_LIBDIR

find_package(mc_rtc REQUIRED)
find_package(qp_solver_collection REQUIRED)
find_package(force_control_collection REQUIRED)
find_package(nmpc_ddp REQUIRED)

if(USE_ROS2)
  # mc_rtc
  find_package(ament_cmake REQUIRED)
  find_package(rclcpp REQUIRED)
else()
  option(BUILD_TESTING "Build test" ON)
endif()

add_subdirectory(src)

if(NOT USE_ROS2)
  install(EXPORT ${PROJECT_NAME}
    FILE ${PROJECT_NAME}Targets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    NAMESPACE ${PROJECT_NAME}::
    )

  include(CMakePackageConfigHelpers)
  # generate the config file that is includes the exports
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

  # generate the version file for the config file
  # Extract version numbers from package.xml
  file(READ package.xml PACKAGE_XML)
  string(REGEX MATCH "<version>[0-9]+\\.[0-9]+\\.[0-9]+</version>" DIRTY_VERSION_STRING ${PACKAGE_XML})
  string(REGEX REPLACE "^<version>([0-9]+)\\.([0-9]+)\\.([0-9]+)</version>$" "\\1" CCC_MAJOR_VERSION "${DIRTY_VERSION_STRING}")
  string(REGEX REPLACE "^<version>([0-9]+)\\.([0-9]+)\\.([0-9]+)</version>$" "\\2" CCC_MINOR_VERSION "${DIRTY_VERSION_STRING}")
  string(REGEX REPLACE "^<version>([0-9]+)\\.([0-9]+)\\.([0-9]+)</version>$" "\\3" CCC_PATCH_VERSION "${DIRTY_VERSION_STRING}")
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${CCC_MAJOR_VERSION}.${CCC_MINOR_VERSION}.${CCC_PATCH_VERSION}"
    COMPATIBILITY AnyNewerVersion
    )

  # install config files
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    )
endif()

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
  ament_export_targets(${PROJECT_NAME})
  ament_package()
endif() 