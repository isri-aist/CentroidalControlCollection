add_library(CCC
  CommonModels.cpp
  PreviewControlZmp.cpp
  DdpZmp.cpp
  DcmTracking.cpp
  FootGuidedControl.cpp
  LinearMpcZmp.cpp
  IntrinsicallyStableMpc.cpp
  SingularPreviewControlZmp.cpp
  LinearMpcZ.cpp
  LinearMpcXY.cpp
  DdpCentroidal.cpp
  )
if(DEFINED CATKIN_DEVEL_PREFIX)
  target_include_directories(CCC PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${catkin_INCLUDE_DIRS}
    )
  target_link_libraries(CCC PUBLIC
    ${catkin_LIBRARIES}
    )
else()
  target_link_libraries(CCC PUBLIC qp_solver_collection::QpSolverCollection nmpc_ddp::nmpc_ddp)
  target_include_directories(CCC PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
  target_compile_definitions(CCC PUBLIC CCC_STANDALONE)
endif()
target_compile_features(CCC PUBLIC cxx_std_17)

if(BUILD_SHARED_LIBS)
  set_target_properties(CCC PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION})
endif()

install(TARGETS CCC
          EXPORT "${TARGETS_EXPORT_NAME}"
          LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/CCC DESTINATION "${INCLUDE_INSTALL_DIR}")
