/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/PreviewControlZmp.h>

#include "SimModels.h"

TEST(TestPreviewControlZmp, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.01; // [sec]
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]

  // Setup preview control
  CCC::PreviewControlZmp pc(com_height, horizon_duration, horizon_dt);

  std::function<double(double)> ref_zmp_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    int phase_idx = static_cast<int>(std::floor(t) - 1);
    double phase_time = t - std::floor(t);
    std::vector<double> zmp_list = {0.0, 0.1, -0.1, 0.1, 0.2, 0.3, 0.4};
    if(phase_idx <= 0)
    {
      return 0.0;
    }
    else if(phase_idx >= zmp_list.size())
    {
      return zmp_list.back();
    }
    else if(phase_time < 0.2)
    {
      double ratio = phase_time / 0.2;
      return (1 - ratio) * zmp_list[phase_idx - 1] + ratio * zmp_list[phase_idx];
    }
    else
    {
      return zmp_list[std::min(phase_idx, static_cast<int>(zmp_list.size()) - 1)];
    }
  };

  // Setup simulation
  ComZmpSimModel sim_model(com_height);
  sim_model.calcDiscMatrix(sim_dt);

  // Setup dump file
  std::string file_path = "/tmp/TestPreviewControlZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos com_vel planned_zmp ref_zmp" << std::endl;

  // Setup control loop
  Eigen::Vector2d com_pos_vel = Eigen::Vector2d::Zero();
  double planned_zmp = com_pos_vel[0];

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    CCC::PreviewControlZmp::InitialParam initial_param;
    initial_param << com_pos_vel, CCC::constants::g / com_height * (com_pos_vel[0] - planned_zmp);
    planned_zmp = pc.planOnce(ref_zmp_func, initial_param, t, sim_dt);

    // Dump
    double ref_zmp = ref_zmp_func(t);
    ofs << t << " " << com_pos_vel.transpose() << " " << planned_zmp << " " << ref_zmp << std::endl;

    // Check
    EXPECT_LT(std::abs(planned_zmp - ref_zmp), 0.1); // [m]

    // Update
    t += sim_dt;
    Eigen::Vector1d sim_input;
    sim_input << planned_zmp;
    com_pos_vel = sim_model.stateEqDisc(com_pos_vel, sim_input);
  }

  // Final check
  double ref_zmp = ref_zmp_func(t);
  EXPECT_LT(std::abs(planned_zmp - ref_zmp), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel[0] - ref_zmp), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel[1]), 1e-2);

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
