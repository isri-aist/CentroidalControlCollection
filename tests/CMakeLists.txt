set(CCC_gtest_list
  TestStateSpaceModel
  TestVariantSequentialExtension
  TestPreviewControlZmp
  TestDdpZmp
  TestDcmTracking
  TestLinearMpcZ
  TestLinearMpcXY
  )

foreach(NAME IN LISTS CCC_gtest_list)
  catkin_add_gtest(${NAME} src/${NAME}.cpp)
  target_link_libraries(${NAME} CCC)
endforeach()
