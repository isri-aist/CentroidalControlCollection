/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

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
  sva::ForceVecd disturb_impulse_per_mass =
      sva::ForceVecd(Eigen::Vector3d::Zero(), Eigen::Vector3d(0.05, 0.05, 0.0)); // [m/s], [m^2/s]

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
      ref_data.pos.linear() << 0.0, 0.0, 1.0; // [m]
    }
    else if(t < 1.6)
    {
      ref_data.pos.linear() << 0.25, 0.0, 1.2; // [m]
    }
    else
    {
      ref_data.pos.linear() << 0.5, 0.0, 1.0; // [m]
    }
    return ref_data;
  };

  // Setup simulation
  CentroidalSim sim(mass, moment_of_inertia, sim_dt);
  sim.state_.pos.linear().z() = 1.0;

  // Setup dump file
  std::string file_path = "/tmp/TestPreviewControlCentroidal.txt";
  std::ofstream ofs(file_path);
  ofs << "time "
         "planned_pos_rx planned_pos_ry planned_pos_rz planned_pos_tx planned_pos_ty planned_pos_tz "
         "ref_pos_rx ref_pos_ry ref_pos_rz ref_pos_tx ref_pos_ty ref_pos_tz "
         "planned_vel_ax planned_vel_ay planned_vel_az planned_vel_lx planned_vel_ly planned_vel_lz "
         "planned_momentum_nx planned_momentum_ny planned_momentum_nz planned_momentum_fx planned_momentum_fy "
         "planned_momentum_fz "
         "wrench_nx wrench_ny wrench_nz wrench_fx wrench_fy wrench_fz "
         "computation_time"
      << std::endl;

  // Run control loop
  sva::ForceVecd planned_wrench =
      sva::ForceVecd(Eigen::Vector3d::Zero(), Eigen::Vector3d(0.0, 0.0, mass * CCC::constants::g));
  constexpr double end_time = 3.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    CCC::PreviewControlCentroidal::InitialParam initial_param;
    initial_param.pos = sim.state_.pos;
    initial_param.vel = sim.state_.vel;
    initial_param.acc.linear() = planned_wrench.force() / mass - Eigen::Vector3d(0, 0, CCC::constants::g);
    initial_param.acc.angular() = planned_wrench.moment().cwiseQuotient(moment_of_inertia);
    planned_wrench = pc.planOnce(motion_param_func(t), ref_data_func, initial_param, t, sim_dt);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    const auto & ref_data = ref_data_func(t);
    ofs << t << " " << sim.state_.pos.vector().transpose() << " " << ref_data.pos.vector().transpose() << " "
        << sim.state_.vel.vector().transpose() << " " << sim.state_.momentum.vector().transpose() << " "
        << planned_wrench.vector().transpose() << " " << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((sim.state_.pos - ref_data.pos).vector().norm(), 2.0); // [m], [rad]
    EXPECT_LT(sim.state_.vel.vector().norm(), 2.0); // [m/s], [rad/s]

    // Simulate
    t += sim_dt;
    sim.update(CentroidalSim::Input(planned_wrench));

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
  EXPECT_LT((sim.state_.pos - ref_data.pos).vector().norm(), 0.1); // [m], [rad]
  EXPECT_LT(sim.state_.vel.vector().norm(), 0.1); // [m/s], [rad/s]

  Eigen::Map<Eigen::VectorXd> computation_duration_vec(computation_duration_list.data(),
                                                       computation_duration_list.size());
  std::cout << "Computation time per control cycle:\n"
            << "  mean: " << computation_duration_vec.mean() << " [ms], stdev: "
            << std::sqrt((computation_duration_vec.array() - computation_duration_vec.mean()).square().mean())
            << " [ms], max: " << computation_duration_vec.maxCoeff() << " [ms]" << std::endl;

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:5 w lp, \"\" u 1:11 w l lw 2 # pos_tx\n"
            << "  plot \"" << file_path << "\" u 1:7 w lp, \"\" u 1:13 w l lw 2 # pos_tz\n"
            << "  plot \"" << file_path << "\" u 1:3 w lp, \"\" u 1:9 w l lw 2 # pos_ry\n"
            << "  plot \"" << file_path << "\" u 1:31 w lp # wrench_fz\n"
            << "  plot \"" << file_path << "\" u 1:32 w lp # computation_time\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
