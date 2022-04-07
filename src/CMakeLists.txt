add_library(CCC
  PreviewControlZmp.cpp
  DdpZmp.cpp
  LinearMpcZ.cpp
  LinearMpcXY.cpp
  )
target_include_directories(CCC PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${catkin_INCLUDE_DIRS}
  )
target_link_libraries(CCC PUBLIC
  ${catkin_LIBRARIES}
  )
