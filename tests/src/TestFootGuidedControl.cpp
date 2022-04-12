/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/FootGuidedControl.h>

#include "FootstepManager.h"
#include "SimModels.h"

TEST(TestFootGuidedControl, Test1)
{
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]

  // Setup foot-guided control
  CCC::FootGuidedControl foot_guided(com_height);

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
  std::string file_path = "/tmp/TestFootGuidedControl.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos_x com_pos_y planned_zmp_x planned_zmp_y ref_zmp_x ref_zmp_y capture_point_x capture_point_y"
      << std::endl;

  // Setup control loop
  Eigen::Vector2d planned_zmp = Eigen::Vector2d::Zero();

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    footstep_manager.update(t);
    CCC::FootGuidedControl::InitialParam initial_param =
        sim.state_.pos() + std::sqrt(com_height / CCC::constants::g) * sim.state_.vel();
    const auto & ref_data = footstep_manager.makeFootGuidedControlRefData(t);
    planned_zmp = foot_guided.planOnce(ref_data, initial_param, t);

    // Dump
    Eigen::Vector2d ref_zmp = footstep_manager.refZmp(t);
    ofs << t << " " << sim.state_.pos().transpose() << " " << planned_zmp.transpose() << " " << ref_zmp.transpose()
        << " " << initial_param.transpose() << std::endl;

    // Check
    EXPECT_LT((planned_zmp - ref_zmp).norm(), 0.1); // [m]

    // Simulate
    t += sim_dt;
    sim.update(planned_zmp);
  }

  // Final check
  Eigen::Vector2d ref_zmp = footstep_manager.refZmp(t);
  EXPECT_LT((planned_zmp - ref_zmp).norm(), 1e-2);
  EXPECT_LT((sim.state_.pos() - ref_zmp).norm(), 1e-2);
  EXPECT_LT(sim.state_.vel().norm(), 1e-2);

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:6 w l lw 2, \"\" u 1:8 w lp # x\n"
            << "  plot \"" << file_path << "\" u 1:3 w lp, \"\" u 1:5 w lp, \"\" u 1:7 w l lw 2, \"\" u 1:9 w lp # y\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
