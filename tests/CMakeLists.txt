set(CCC_gtest_list
  testStateSpaceModel
  testVariantSequentialExtension
  )

foreach(NAME IN LISTS CCC_gtest_list)
  catkin_add_gtest(${NAME} src/${NAME}.cpp)
  target_link_libraries(${NAME} CCC)
endforeach()
