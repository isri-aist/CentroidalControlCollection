if(NOT DEFINED CATKIN_DEVEL_PREFIX)
  find_package(GTest REQUIRED)
  include(GoogleTest)
  function(add_CCC_test NAME)
    add_executable(${NAME} src/${NAME}.cpp)
    target_link_libraries(${NAME} PUBLIC GTest::gtest CCC)
    gtest_discover_tests(${NAME})
  endfunction()
else()
  function(add_CCC_test NAME)
    catkin_add_gtest(${NAME} src/${NAME}.cpp)
    target_link_libraries(${NAME} CCC)
  endfunction()
endif()

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
  TestLinearMpcZ
  TestLinearMpcXY
  TestPreviewControlCentroidal
  TestDdpCentroidal
  )

foreach(NAME IN LISTS CCC_gtest_list)
  add_CCC_test(${NAME})
endforeach()
