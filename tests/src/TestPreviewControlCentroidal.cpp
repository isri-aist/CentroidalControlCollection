/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <ForceColl/Contact.h>

#include <CCC/Constants.h>
#include <CCC/PreviewControlCentroidal.h>

#include "ContactManager.h"
#include "SimModels.h"

TEST(TestPreviewControlCentroidal, PlanOnce)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.01; // [sec]
  double sim_dt = 0.005; // [sec]
  double mass = 100.0; // [kg]
  Eigen::Vector3d moment_of_inertia = Eigen::Vector3d(40.0, 20.0, 10.0); // [kg m^2]
  std::vector<double> disturb_time_list = {1.0}; // [sec]
  Eigen::Vector6d disturb_impulse_per_mass =
      (Eigen::Vector6d() << 0.05, 0.05, 0.0, 0.0, 0.0, 0.0).finished(); // [m/s], [m^2/s]

  // Setup preview control
  std::vector<double> computation_duration_list;
  CCC::PreviewControlCentroidal pc(mass, moment_of_inertia, horizon_duration, horizon_dt);

  // Setup contact
  std::function<CCC::PreviewControlCentroidal::MotionParam(double)> motion_param_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::PreviewControlCentroidal::MotionParam motion_param;
    if(t < 1.4)
    {
      motion_param.contact_list.push_back(
          makeContactFromRect({Eigen::Vector2d(-0.1, -0.1), Eigen::Vector2d(0.1, 0.1)}));
    }
    else if(t < 1.6)
    {
      motion_param.contact_list.push_back(
          makeContactFromRect({Eigen::Vector2d(0.15, 0.15), Eigen::Vector2d(0.35, 0.35)}));
    }
    else
    {
      motion_param.contact_list.push_back(makeContactFromRect({Eigen::Vector2d(0.4, -0.1), Eigen::Vector2d(0.6, 0.1)}));
    }
    return motion_param;
  };
  std::function<CCC::PreviewControlCentroidal::RefData(double)> ref_data_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::PreviewControlCentroidal::RefData ref_data;
    if(t < 1.4)
    {
      ref_data.pos.head<3>() << 0.0, 0.0, 1.0; // [m]
    }
    else if(t < 1.6)
    {
      ref_data.pos.head<3>() << 0.25, 0.0, 1.2; // [m]
    }
    else
    {
      ref_data.pos.head<3>() << 0.5, 0.0, 1.0; // [m]
    }
    return ref_data;
  };

  // Setup simulation
  CentroidalSim sim(mass, moment_of_inertia, sim_dt);
  sim.state_.pos.head<3>().z() = 1.0;

  // Setup dump file
  std::string file_path = "/tmp/TestPreviewControlCentroidal.txt";
  std::ofstream ofs(file_path);
  ofs << "time planned_com_pos_x planned_com_pos_y planned_com_pos_z planned_ori_x planned_ori_y planned_ori_z "
         "ref_com_pos_x ref_com_pos_y ref_com_pos_z ref_ori_x ref_ori_y ref_ori_z planned_com_vel_x planned_com_vel_y "
         "planned_com_vel_z planned_angular_vel_x planned_angular_vel_y planned_angular_vel_z "
         "planned_linear_momentum_x planned_linear_momentum_y planned_linear_momentum_z planned_angular_momentum_x "
         "planned_angular_momentum_y planned_angular_momentum_z wrench_fx wrench_fy wrench_fz wrench_nx wrench_ny "
         "wrench_nz computation_time"
      << std::endl;

  // Run control loop
  Eigen::Vector6d planned_wrench = (Eigen::Vector6d() << 0.0, 0.0, mass * CCC::constants::g, 0.0, 0.0, 0.0).finished();
  constexpr double end_time = 3.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    CCC::PreviewControlCentroidal::InitialParam initial_param;
    initial_param.pos = sim.state_.pos;
    initial_param.vel = sim.state_.vel;
    initial_param.acc.head<3>() = planned_wrench.head<3>() / mass - Eigen::Vector3d(0, 0, CCC::constants::g);
    initial_param.acc.tail<3>() = planned_wrench.tail<3>().cwiseQuotient(moment_of_inertia);
    planned_wrench = pc.planOnce(motion_param_func(t), ref_data_func, initial_param, t, sim_dt);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    const auto & ref_data = ref_data_func(t);
    ofs << t << " " << sim.state_.pos.transpose() << " " << ref_data.pos.transpose() << " "
        << sim.state_.vel.transpose() << " " << sim.state_.momentum.transpose() << " " << planned_wrench.transpose()
        << " " << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((sim.state_.pos - ref_data.pos).norm(), 2.0); // [m], [rad]
    EXPECT_LT(sim.state_.vel.norm(), 2.0); // [m/s], [rad/s]

    // Simulate
    t += sim_dt;
    sim.update(planned_wrench);

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
  const auto & ref_data = ref_data_func(t);
  EXPECT_LT((sim.state_.pos - ref_data.pos).norm(), 0.1); // [m], [rad]
  EXPECT_LT(sim.state_.vel.norm(), 0.1); // [m/s], [rad/s]

  Eigen::Map<Eigen::VectorXd> computation_duration_vec(computation_duration_list.data(),
                                                       computation_duration_list.size());
  std::cout << "Computation time per control cycle:\n"
            << "  mean: " << computation_duration_vec.mean() << " [ms], stdev: "
            << std::sqrt((computation_duration_vec.array() - computation_duration_vec.mean()).square().mean())
            << " [ms], max: " << computation_duration_vec.maxCoeff() << " [ms]" << std::endl;

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:8 w l lw 2 # pos_x\n"
            << "  plot \"" << file_path << "\" u 1:4 w lp, \"\" u 1:10 w l lw 2 # pos_z\n"
            << "  plot \"" << file_path << "\" u 1:6 w lp, \"\" u 1:12 w l lw 2 # ori_y\n"
            << "  plot \"" << file_path << "\" u 1:28 w lp # force_z\n"
            << "  plot \"" << file_path << "\" u 1:32 w lp # computation_time\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
