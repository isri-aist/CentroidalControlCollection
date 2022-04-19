/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DcmTracking.h>

#include "FootstepManager.h"
#include "SimModels.h"

TEST(TestDcmTracking, Test1)
{
  double sim_dt = 0.005; // [sec]
  double com_height = 1.0; // [m]
  std::vector<double> disturb_time_list = {4.5, 8.5}; // [sec]
  Eigen::Vector2d disturb_impulse_per_mass = Eigen::Vector2d(0.05, 0.05); // [m/s]

  // Setup DCM tracking
  std::vector<double> computation_duration_list;
  CCC::DcmTracking dcm_tracking(com_height);

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
  std::string file_path = "/tmp/TestDcmTracking.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos_x com_pos_y planned_zmp_x planned_zmp_y ref_zmp_x ref_zmp_y capture_point_x capture_point_y "
         "computation_time"
      << std::endl;

  // Setup control loop
  Eigen::Vector2d planned_zmp = Eigen::Vector2d::Zero();

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    footstep_manager.update(t);
    CCC::DcmTracking::InitialParam initial_param =
        sim.state_.pos() + std::sqrt(com_height / CCC::constants::g) * sim.state_.vel();
    const auto & ref_data = footstep_manager.makeDcmTrackingRefData(t);
    planned_zmp = dcm_tracking.planOnce(ref_data, initial_param, t);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    Eigen::Vector2d ref_zmp = ref_data.current_zmp;
    ofs << t << " " << sim.state_.pos().transpose() << " " << planned_zmp.transpose() << " " << ref_zmp.transpose()
        << " " << initial_param.transpose() << " " << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((planned_zmp - ref_zmp).norm(), 0.1); // [m]

    // Simulate
    t += sim_dt;
    sim.update(planned_zmp);

    // Add disturbance
    for(double disturb_time : disturb_time_list)
    {
      if(disturb_time <= t && t < disturb_time + sim_dt)
      {
        sim.addDisturb(disturb_impulse_per_mass);
        break;
      }
    }
  }

  // Final check
  Eigen::Vector2d ref_zmp = footstep_manager.makeDcmTrackingRefData(t).current_zmp;
  EXPECT_LT((planned_zmp - ref_zmp).norm(), 1e-2);
  EXPECT_LT((sim.state_.pos() - ref_zmp).norm(), 1e-2);
  EXPECT_LT(sim.state_.vel().norm(), 1e-2);

  Eigen::Map<Eigen::VectorXd> computation_duration_vec(computation_duration_list.data(),
                                                       computation_duration_list.size());
  std::cout << "Computation time per control cycle:\n"
            << "  mean: " << computation_duration_vec.mean() << " [ms], stdev: "
            << std::sqrt((computation_duration_vec.array() - computation_duration_vec.mean()).square().mean())
            << " [ms], max: " << computation_duration_vec.maxCoeff() << " [ms]" << std::endl;

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:6 w l lw 2 # x\n"
            << "  plot \"" << file_path << "\" u 1:3 w lp, \"\" u 1:5 w lp, \"\" u 1:7 w l lw 2 # y\n"
            << "  plot \"" << file_path << "\" u 1:10 w lp # computation_time\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
