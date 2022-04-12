/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/LinearMpcZmp.h>

#include "FootstepManager.h"
#include "SimModels.h"

TEST(TestLinearMpcZmp, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.02; // [sec]
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]

  // Setup MPC
  CCC::LinearMpcZmp mpc(com_height, horizon_duration, horizon_dt);

  // Setup footstep
  FootstepManager footstep_manager;
  double transit_duration = 0.2; // [sec]
  double swing_duration = 0.8; // [sec]
  footstep_manager.appendFootstep(
      Footstep(Foot::Left, Eigen::Vector2d(0.2, 0.1), 2.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Right, Eigen::Vector2d(0.4, -0.1), 3.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Left, Eigen::Vector2d(0.6, 0.1), 4.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Right, Eigen::Vector2d(0.8, -0.1), 5.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Left, Eigen::Vector2d(0.6, 0.1), 6.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Right, Eigen::Vector2d(0.6, -0.1), 7.0, transit_duration, swing_duration));

  // Setup simulation
  ComZmpSim2d sim(com_height, sim_dt);

  // Setup dump file
  std::string file_path = "/tmp/TestLinearMpcZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos_x com_pos_y planned_zmp_x planned_zmp_y ref_zmp_min_x ref_zmp_min_y ref_zmp_max_x ref_zmp_max_y"
      << std::endl;

  // Setup control loop
  Eigen::Vector2d planned_zmp = sim.state_.pos();

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    footstep_manager.update(t);
    CCC::LinearMpcZmp::InitialParam initial_param;
    initial_param.pos = sim.state_.pos();
    initial_param.vel = sim.state_.vel();
    initial_param.acc = CCC::constants::g / com_height * (sim.state_.pos() - planned_zmp);
    planned_zmp =
        mpc.planOnce(std::bind(&FootstepManager::makeLinearMpcZmpRefData, &footstep_manager, std::placeholders::_1),
                     initial_param, t, sim_dt);

    // Dump
    const auto & ref_data = footstep_manager.makeLinearMpcZmpRefData(t);
    ofs << t << " " << sim.state_.pos().transpose() << " " << planned_zmp.transpose() << " "
        << ref_data.zmp_limits[0].transpose() << " " << ref_data.zmp_limits[1].transpose() << std::endl;

    // Check
    EXPECT_TRUE(((planned_zmp - ref_data.zmp_limits[0]).array() >= 0).all());
    EXPECT_TRUE(((ref_data.zmp_limits[1] - planned_zmp).array() >= 0).all());

    // Simulate
    t += sim_dt;
    sim.update(planned_zmp);
  }

  // Final check
  const auto & ref_data = footstep_manager.makeLinearMpcZmpRefData(t);
  EXPECT_TRUE(((planned_zmp - ref_data.zmp_limits[0]).array() >= 0).all());
  EXPECT_TRUE(((ref_data.zmp_limits[1] - planned_zmp).array() >= 0).all());
  EXPECT_TRUE(((sim.state_.pos() - ref_data.zmp_limits[0]).array() >= 0).all());
  EXPECT_TRUE(((ref_data.zmp_limits[1] - sim.state_.pos()).array() >= 0).all());

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path
            << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:6 w l lw 2, \"\" u 1:8 w l lw 2 # x\n"
            << "  plot \"" << file_path
            << "\" u 1:3 w lp, \"\" u 1:5 w lp, \"\" u 1:7 w l lw 2, \"\" u 1:9 w l lw 2 # y\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
