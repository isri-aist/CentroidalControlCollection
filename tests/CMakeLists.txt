find_package(Boost REQUIRED COMPONENTS unit_test_framework)

set(ccc_TEST
  testStateSpaceModel
)

foreach(NAME IN LISTS ccc_TEST)
  add_executable(${NAME} src/${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC Boost::unit_test_framework Boost::disable_autolinking CCC)
  target_compile_definitions(${NAME} PRIVATE -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)
  target_include_directories(${NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
  add_test(${NAME} ${NAME})
endforeach()
