find_package(GTest REQUIRED)
include(GoogleTest)
# Prevents discovery failure before install
# No effect in CMake < 3.18
set(CMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE PRE_TEST)

function(add_CCC_test NAME)
  add_executable(${NAME} src/${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC GTest::gtest CCC)
  gtest_discover_tests(${NAME})
endfunction()

set(CCC_gtest_list
  TestStateSpaceModel
  TestVariantSequentialExtension
  TestInvariantSequentialExtension
  TestPreviewControlZmp
  TestDdpZmp
  TestDcmTracking
  TestFootGuidedControl
  TestLinearMpcZmp
  TestIntrinsicallyStableMpc
  TestSingularPreviewControlZmp
  TestStepMpc
  TestLinearMpcZ
  TestLinearMpcXY
  TestPreviewControlCentroidal
  TestDdpCentroidal
  TestDdpSingleRigidBody
  )

foreach(NAME IN LISTS CCC_gtest_list)
  add_CCC_test(${NAME})
endforeach()

if(ENABLE_PYBULLET_TEST)
  if(USE_ROS2)
    find_package(launch_testing_ament_cmake)
    ament_add_gtest(TestSimDdpSingleRigidBody src/TestSimDdpSingleRigidBody.cpp TIMEOUT 600)
    ament_target_dependencies(TestSimDdpSingleRigidBody
      rclcpp
      std_msgs
      std_srvs
    )

    target_link_libraries(TestSimDdpSingleRigidBody CCC)
    install(TARGETS TestSimDdpSingleRigidBody DESTINATION lib/${PROJECT_NAME})
    add_launch_test(scripts/simSingleRigidBody.py)
    install(DIRECTORY scripts DESTINATION share/${PROJECT_NAME}/tests)
  endif()
endif()