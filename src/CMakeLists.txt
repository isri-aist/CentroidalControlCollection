add_library(CCC
  LinearMpcZ.cpp
  LinearMpcXY.cpp
  )
target_include_directories(CCC PUBLIC
  ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(CCC PUBLIC
  ${catkin_LIBRARIES}
  )
