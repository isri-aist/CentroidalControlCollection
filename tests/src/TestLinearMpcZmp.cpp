/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/LinearMpcZmp.h>

#include "SimModels.h"

TEST(TestLinearMpcZmp, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.04; // [sec]
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]

  // Setup MPC
  CCC::LinearMpcZmp mpc(com_height, horizon_duration, horizon_dt);

  std::function<CCC::LinearMpcZmp::RefData(double)> ref_data_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::LinearMpcZmp::RefData ref_data;
    int phase_idx = static_cast<int>(std::floor(t) - 1);
    double phase_time = t - std::floor(t);
    std::vector<double> zmp_list = {0.0, 0.1, -0.1, 0.1, 0.2, 0.3, 0.4};
    double support_region_width = 0.05; // [m]
    if(phase_idx <= 0)
    {
      double ref_zmp = 0.0;
      ref_data.zmp_limits = {ref_zmp - 0.5 * support_region_width, ref_zmp + 0.5 * support_region_width};
    }
    else if(phase_idx >= zmp_list.size())
    {
      double ref_zmp = zmp_list.back();
      ref_data.zmp_limits = {ref_zmp - 0.5 * support_region_width, ref_zmp + 0.5 * support_region_width};
    }
    else if(phase_time < 0.2)
    {
      double ref_zmp_min = std::min(zmp_list[phase_idx - 1], zmp_list[phase_idx]);
      double ref_zmp_max = std::max(zmp_list[phase_idx - 1], zmp_list[phase_idx]);
      ref_data.zmp_limits = {ref_zmp_min - 0.5 * support_region_width, ref_zmp_max + 0.5 * support_region_width};
    }
    else
    {
      double ref_zmp = zmp_list[std::min(phase_idx, static_cast<int>(zmp_list.size()) - 1)];
      ref_data.zmp_limits = {ref_zmp - 0.5 * support_region_width, ref_zmp + 0.5 * support_region_width};
    }
    return ref_data;
  };

  // Setup simulation
  ComZmpSimModel sim_model(com_height);
  sim_model.calcDiscMatrix(sim_dt);

  // Setup dump file
  std::string file_path = "/tmp/TestLinearMpcZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos com_vel planned_zmp ref_zmp_min ref_zmp_max" << std::endl;

  // Setup control loop
  Eigen::Vector2d com_pos_vel = Eigen::Vector2d::Zero();
  double planned_zmp = com_pos_vel[0];

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    CCC::LinearMpcZmp::InitialParam initial_param;
    initial_param << com_pos_vel, CCC::constants::g / com_height * (com_pos_vel[0] - planned_zmp);
    planned_zmp = mpc.planOnce(ref_data_func, initial_param, t, sim_dt);

    // Dump
    const auto & ref_data = ref_data_func(t);
    ofs << t << " " << com_pos_vel.transpose() << " " << planned_zmp << " " << ref_data.zmp_limits[0] << " "
        << ref_data.zmp_limits[1] << std::endl;

    // Check
    EXPECT_LE(ref_data.zmp_limits[0], planned_zmp);
    EXPECT_LE(planned_zmp, ref_data.zmp_limits[1]);

    // Simulate
    t += sim_dt;
    Eigen::Vector1d sim_input;
    sim_input << planned_zmp;
    com_pos_vel = sim_model.stateEqDisc(com_pos_vel, sim_input);
  }

  // Final check
  const auto & ref_data = ref_data_func(t);
  EXPECT_LE(ref_data.zmp_limits[0], planned_zmp);
  EXPECT_LE(planned_zmp, ref_data.zmp_limits[1]);
  EXPECT_LE(ref_data.zmp_limits[0], com_pos_vel[0]);
  EXPECT_LE(com_pos_vel[0], ref_data.zmp_limits[1]);

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:5 w l lw 2, \"\" u 1:6 w l lw 2\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
