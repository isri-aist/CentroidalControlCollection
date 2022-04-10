/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DcmTracking.h>
#include <CCC/EigenTypes.h>

#include "SimModels.h"

/** \brief CoM position and velocity. */
struct ComPosVel
{
  Eigen::Vector2d x = Eigen::Vector2d::Zero();
  Eigen::Vector2d y = Eigen::Vector2d::Zero();
};

TEST(TestDcmTracking, Test1)
{
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]
  double omega = std::sqrt(CCC::constants::g / com_height);

  // Setup DCM tracking
  CCC::DcmTracking dcm_tracking(com_height);

  // clang-format off
  std::map<double, Eigen::Vector2d> time_zmp_list_all = {
      {0.0, Eigen::Vector2d::Constant(0.0)},
      {2.0, Eigen::Vector2d::Constant(0.1)},
      {3.0, Eigen::Vector2d::Constant(-0.1)},
      {4.0, Eigen::Vector2d::Constant(0.1)},
      {5.0, Eigen::Vector2d::Constant(0.2)},
      {6.0, Eigen::Vector2d::Constant(0.3)},
      {7.0, Eigen::Vector2d::Constant(0.4)}
  };
  // clang-format on

  // Setup simulation
  ComZmpSimModel sim_model(com_height);
  sim_model.calcDiscMatrix(sim_dt);

  // Setup dump file
  std::string file_path = "/tmp/TestDcmTracking.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos com_vel planned_zmp ref_zmp" << std::endl;

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
    CCC::DcmTracking::InitialParam initial_param;
    initial_param << com_pos_vel.x[0] + com_pos_vel.x[1] / omega, com_pos_vel.y[0] + com_pos_vel.y[1] / omega;
    CCC::DcmTracking::RefData ref_data;
    ref_data.current_zmp = std::prev(time_zmp_list_all.upper_bound(t))->second;
    std::map<double, Eigen::Vector2d> time_zmp_list_future;
    constexpr int time_zmp_num = 3;
    for(auto time_zmp_it = time_zmp_list_all.upper_bound(t);
        time_zmp_it != time_zmp_list_all.end() && time_zmp_list_future.size() < time_zmp_num; time_zmp_it++)
    {
      time_zmp_list_future.emplace(*time_zmp_it);
    }
    ref_data.time_zmp_list = time_zmp_list_future;
    planned_zmp = dcm_tracking.planOnce(ref_data, initial_param, t);

    // Dump
    ref_zmp = ref_data.current_zmp;
    ofs << t << " " << com_pos_vel.x.transpose() << " " << planned_zmp.x() << " " << ref_zmp.x() << std::endl;

    // Check
    EXPECT_LT(std::abs(planned_zmp.x() - ref_zmp.x()), 0.1); // [m]

    // Simulate
    t += sim_dt;
    Eigen::Vector1d sim_input;
    sim_input << planned_zmp.x();
    com_pos_vel.x = sim_model.stateEqDisc(com_pos_vel.x, sim_input);
    sim_input << planned_zmp.y();
    com_pos_vel.y = sim_model.stateEqDisc(com_pos_vel.y, sim_input);
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
