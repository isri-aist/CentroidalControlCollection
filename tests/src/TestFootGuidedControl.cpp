/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/FootGuidedControl.h>

#include "SimModels.h"

/** \brief CoM position and velocity. */
struct ComPosVel
{
  Eigen::Vector2d x = Eigen::Vector2d::Zero();
  Eigen::Vector2d y = Eigen::Vector2d::Zero();
};

TEST(TestFootGuidedControl, Test1)
{
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]
  double omega = std::sqrt(CCC::constants::g / com_height);

  // Setup DCM tracking
  CCC::FootGuidedControl foot_guided(com_height);

  // clang-format off
  std::map<double, Eigen::Vector2d> time_zmp_list_all = {
      {0.0, Eigen::Vector2d::Constant(0.0)},
      {2.0, Eigen::Vector2d::Constant(0.1)},
      {3.5, Eigen::Vector2d::Constant(-0.1)},
      {5.0, Eigen::Vector2d::Constant(0.1)},
      {6.5, Eigen::Vector2d::Constant(0.2)},
      {8.0, Eigen::Vector2d::Constant(0.3)}
  };
  // clang-format on

  // Setup simulation
  ComZmpSimModel sim_model(com_height);
  sim_model.calcDiscMatrix(sim_dt);

  // Setup dump file
  std::string file_path = "/tmp/TestFootGuidedControl.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos capture_point planned_zmp ref_zmp" << std::endl;

  // Setup control loop
  ComPosVel com_pos_vel;
  Eigen::Vector2d ref_zmp = Eigen::Vector2d::Zero();
  Eigen::Vector2d planned_zmp = Eigen::Vector2d::Zero();

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    CCC::FootGuidedControl::InitialParam initial_param;
    initial_param = com_pos_vel.x[0] + com_pos_vel.x[1] / omega;
    CCC::FootGuidedControl::RefData ref_data;
    auto time_zmp_it = time_zmp_list_all.upper_bound(t + sim_dt);
    ref_data.transit_start_zmp = std::prev(time_zmp_it)->second.x();
    if(time_zmp_it == time_zmp_list_all.end())
    {
      double transit_offset_duration = 1.0; // [sec]
      ref_data.transit_end_zmp = ref_data.transit_start_zmp;
      ref_data.transit_start_time = t + transit_offset_duration;
      ref_zmp.x() = ref_data.transit_start_zmp;
    }
    else
    {
      ref_data.transit_end_zmp = time_zmp_it->second.x();
      double transit_duration = 0.2; // [sec]
      ref_data.transit_start_time = time_zmp_it->first - transit_duration;
      ref_data.transit_duration = transit_duration;
      if(t < ref_data.transit_start_time)
      {
        ref_zmp.x() = ref_data.transit_start_zmp;
      }
      else
      {
        double zmp_transit_vel = (ref_data.transit_end_zmp - ref_data.transit_start_zmp) / transit_duration;
        ref_zmp.x() = ref_data.transit_start_zmp + zmp_transit_vel * (t - ref_data.transit_start_time);
      }
    }
    planned_zmp.x() = foot_guided.planOnce(ref_data, initial_param, t);

    // Dump
    ofs << t << " " << com_pos_vel.x[0] << " " << initial_param << " " << planned_zmp.x() << " " << ref_zmp.x()
        << std::endl;

    // Check
    EXPECT_LT(std::abs(planned_zmp.x() - ref_zmp.x()), 0.1); // [m]

    // Simulate
    t += sim_dt;
    Eigen::Vector1d sim_input;
    sim_input << planned_zmp.x();
    com_pos_vel.x = sim_model.stateEqDisc(com_pos_vel.x, sim_input);
  }

  // Final check
  EXPECT_LT(std::abs(planned_zmp.x() - ref_zmp.x()), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel.x[0] - ref_zmp.x()), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel.x[1]), 1e-2);

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:5 w l lw 2\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
