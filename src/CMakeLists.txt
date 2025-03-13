add_library(CCC
  CommonModels.cpp
  PreviewControlZmp.cpp
  DdpZmp.cpp
  DcmTracking.cpp
  FootGuidedControl.cpp
  LinearMpcZmp.cpp
  IntrinsicallyStableMpc.cpp
  SingularPreviewControlZmp.cpp
  StepMpc.cpp
  LinearMpcZ.cpp
  LinearMpcXY.cpp
  PreviewControlCentroidal.cpp
  DdpCentroidal.cpp
  DdpSingleRigidBody.cpp
)

target_include_directories(CCC PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

if(USE_ROS2)
  target_link_libraries(CCC PUBLIC rclcpp::rclcpp)
else()
  target_compile_definitions(CCC PUBLIC CCC_STANDALONE)
endif()

target_link_libraries(CCC PUBLIC
  mc_rtc::mc_rtc_utils 
  mc_rtc::mc_rtc_gui
  qp_solver_collection::QpSolverCollection
  force_control_collection::ForceColl
  nmpc_ddp::nmpc_ddp
)

target_compile_features(CCC PUBLIC cxx_std_17)

if(BUILD_SHARED_LIBS)
  set_target_properties(CCC PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION})
endif()

install(TARGETS CCC
  EXPORT CCC
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/CCC DESTINATION include)
