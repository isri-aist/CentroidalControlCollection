/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/LinearMpcXY.h>

#include "ContactManager.h"
#include "SimModels.h"

TEST(TestLinearMpcXY, Test1)
{
  double horizon_duration = 1.5; // [sec]
  double horizon_dt = 0.1; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.05; // [sec]
  double mass = 100.0; // [kg]
  Eigen::Vector3d moment_of_inertia = Eigen::Vector3d(40.0, 20.0, 10.0); // [kg m^2]

  // Setup MPC
  std::vector<double> computation_duration_list;
  CCC::LinearMpcXY mpc(mass, horizon_dt, horizon_steps);

  // Setup contact
  std::function<CCC::LinearMpcXY::MotionParam(double)> motion_param_func = [mass](double t)
  {
    CCC::LinearMpcXY::MotionParam motion_param;
    motion_param.com_z = 1.0; // [m]
    motion_param.total_force_z = mass * CCC::constants::g; // [N]
    std::array<Eigen::Vector2d, 2> rect_min_max;
    if(t < 3.0)
    {
      rect_min_max = {Eigen::Vector2d(0.9, -0.15), Eigen::Vector2d(1.1, 0.15)}; // [m]
    }
    else if(t < 4.0)
    {
      rect_min_max = {Eigen::Vector2d(0.9, 0.05), Eigen::Vector2d(1.1, 0.15)};
    }
    else if(t < 5.0)
    {
      rect_min_max = {Eigen::Vector2d(1.15, -0.15), Eigen::Vector2d(1.35, -0.05)};
    }
    else if(t < 6.0)
    {
      rect_min_max = {Eigen::Vector2d(1.4, 0.05), Eigen::Vector2d(1.6, 0.15)};
    }
    else
    {
      rect_min_max = {Eigen::Vector2d(1.4, -0.15), Eigen::Vector2d(1.6, 0.15)};
    }
    motion_param.contact_list.push_back(makeContactFromRect(rect_min_max));
    return motion_param;
  };
  std::function<CCC::LinearMpcXY::RefData(double)> ref_data_func = [](double t)
  {
    CCC::LinearMpcXY::RefData ref_data;
    if(t < 3.0)
    {
      ref_data.pos << 1.0, 0.0; // [m]
    }
    else if(t < 4.0)
    {
      ref_data.pos << 1.0, 0.1;
    }
    else if(t < 5.0)
    {
      ref_data.pos << 1.25, -0.1;
    }
    else if(t < 6.0)
    {
      ref_data.pos << 1.5, 0.1;
    }
    else
    {
      ref_data.pos << 1.5, 0.0;
    }
    return ref_data;
  };

  // Setup simulation
  CentroidalSim sim(mass, moment_of_inertia, sim_dt);
  sim.state_.pos.linear() << ref_data_func(0.0).pos, motion_param_func(0.0).com_z;

  // Setup dump file
  std::string file_path = "/tmp/TestLinearMpcXY.txt";
  std::ofstream ofs(file_path);
  ofs << "time planned_com_pos_x planned_com_pos_y planned_com_pos_z ref_com_pos_x ref_com_pos_y ref_com_pos_z "
         "planned_com_vel_x planned_com_vel_y planned_com_vel_z planned_angular_momentum_x planned_angular_momentum_y "
         "planned_angular_momentum_z wrench_nx wrench_ny wrench_nz wrench_fx wrench_fy wrench_fz "
         "computation_time"
      << std::endl;

  // Run control loop
  constexpr double end_time = 8.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    CCC::LinearMpcXY::InitialParam initial_param;
    initial_param.pos = sim.state_.pos.linear().head<2>();
    initial_param.vel = sim.state_.vel.linear().head<2>();
    initial_param.angular_momentum = sim.state_.momentum.moment().head<2>();
    Eigen::VectorXd planned_force_scales = mpc.planOnce(motion_param_func, ref_data_func, initial_param, t);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    const auto & motion_param = motion_param_func(t);
    const auto & ref_data = ref_data_func(t);
    Eigen::Vector3d ref_pos;
    ref_pos << ref_data.pos, motion_param.com_z;
    const auto & planned_total_wrench =
        ForceColl::calcTotalWrench(motion_param.contact_list, planned_force_scales, sim.state_.pos.linear());
    ofs << t << " " << sim.state_.pos.linear().transpose() << " " << ref_pos.transpose() << " "
        << sim.state_.vel.linear().transpose() << " " << sim.state_.momentum.moment().transpose() << " "
        << planned_total_wrench.vector().transpose() << " " << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((sim.state_.pos.linear() - ref_pos).norm(), 2.0); // [m]
    EXPECT_LT(sim.state_.vel.linear().norm(), 2.0); // [m/s]
    EXPECT_LT(sim.state_.momentum.moment().norm(), 5.0); // [kg m^2/s]

    // Simulate
    t += sim_dt;
    sim.update(CentroidalSim::Input(planned_total_wrench));
  }

  // Final check
  const auto & motion_param = motion_param_func(t);
  const auto & ref_data = ref_data_func(t);
  Eigen::Vector3d ref_pos;
  ref_pos << ref_data.pos, motion_param.com_z;
  EXPECT_LT((sim.state_.pos.linear() - ref_pos).norm(), 0.1); // [m]
  EXPECT_LT(sim.state_.vel.linear().norm(), 0.1); // [m/s]
  EXPECT_LT(sim.state_.momentum.moment().norm(), 0.1); // [kg m^2/s]

  Eigen::Map<Eigen::VectorXd> computation_duration_vec(computation_duration_list.data(),
                                                       computation_duration_list.size());
  std::cout << "Computation time per control cycle:\n"
            << "  mean: " << computation_duration_vec.mean() << " [ms], stdev: "
            << std::sqrt((computation_duration_vec.array() - computation_duration_vec.mean()).square().mean())
            << " [ms], max: " << computation_duration_vec.maxCoeff() << " [ms]" << std::endl;

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:5 w l lw 2 # pos_x\n"
            << "  plot \"" << file_path << "\" u 1:3 w lp, \"\" u 1:6 w l lw 2 # pos_y\n"
            << "  plot \"" << file_path << "\" u 1:4 w lp, \"\" u 1:7 w l lw 2 # pos_z\n"
            << "  plot \"" << file_path << "\" u 1:12 w lp # angular_momentum_y\n"
            << "  plot \"" << file_path << "\" u 1:17 w lp # force_x\n"
            << "  plot \"" << file_path << "\" u 1:20 w lp # computation_time\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
