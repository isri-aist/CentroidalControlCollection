add_library(CCC
  CommonModels.cpp
  PreviewControlZmp.cpp
  DdpZmp.cpp
  DcmTracking.cpp
  FootGuidedControl.cpp
  LinearMpcZmp.cpp
  IntrinsicallyStableMpc.cpp
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
