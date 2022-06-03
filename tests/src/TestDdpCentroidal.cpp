/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DdpCentroidal.h>

#include "ContactManager.h"
#include "SimModels.h"

TEST(TestDdpCentroidal, Test1)
{
  double horizon_dt = 0.03; // [sec]
  double horizon_duration = 3.0; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.005; // [sec]
  double mass = 100.0; // [kg]
  std::vector<double> disturb_time_list = {1.0}; // [sec]
  Eigen::Vector3d disturb_impulse_per_mass = Eigen::Vector3d(0.05, 0.05, 0.0); // [m/s]

  // Setup DDP
  std::vector<double> computation_duration_list;
  CCC::DdpCentroidal ddp(mass, horizon_dt, horizon_steps);
  ddp.ddp_solver_->config().max_iter = 3;

  // Setup contact
  std::function<CCC::DdpCentroidal::MotionParam(double)> motion_param_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::DdpCentroidal::MotionParam motion_param;
    if(t < 1.4)
    {
      motion_param.vertex_ridge_list =
          makeVertexRidgeListFromRect({Eigen::Vector2d(-0.1, -0.1), Eigen::Vector2d(0.1, 0.1)});
    }
    else if(t < 1.6)
    {
      motion_param.vertex_ridge_list.setZero(6, 0);
    }
    else
    {
      motion_param.vertex_ridge_list =
          makeVertexRidgeListFromRect({Eigen::Vector2d(0.4, -0.1), Eigen::Vector2d(0.6, 0.1)});
    }
    return motion_param;
  };
  std::function<CCC::DdpCentroidal::RefData(double)> ref_data_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::DdpCentroidal::RefData ref_data;
    if(t < 1.5)
    {
      ref_data.pos << 0.0, 0.0, 1.0; // [m]
    }
    else
    {
      ref_data.pos << 0.5, 0.0, 1.0; // [m]
    }
    return ref_data;
  };

  // Setup simulation
  CentroidalSim sim(mass, sim_dt);
  sim.state_.pos << 0.0, 0.0, 1.0;

  // Setup dump file
  std::string file_path = "/tmp/TestDdpCentroidal.txt";
  std::ofstream ofs(file_path);
  ofs << "time planned_com_pos_x planned_com_pos_y planned_com_pos_z ref_com_pos_x ref_com_pos_y ref_com_pos_z "
         "planned_com_vel_x planned_com_vel_y planned_com_vel_z planned_angular_momentum_x planned_angular_momentum_y "
         "planned_angular_momentum_z wrench_fx wrench_fy wrench_fz wrench_nx wrench_ny wrench_nz ddp_iter "
         "computation_time"
      << std::endl;

  // Run control loop
  constexpr double end_time = 3.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    CCC::DdpCentroidal::InitialParam initial_param;
    initial_param.pos = sim.state_.pos;
    initial_param.vel = sim.state_.vel;
    initial_param.angular_momentum = sim.state_.angular_momentum;
    initial_param.u_list = ddp.ddp_solver_->controlData().u_list;
    if(!initial_param.u_list.empty())
    {
      for(int i = 0; i < ddp.ddp_solver_->config().horizon_steps; i++)
      {
        double tmp_time = t + i * ddp.ddp_problem_->dt();
        int input_dim = ddp.ddp_problem_->inputDim(tmp_time);
        if(initial_param.u_list[i].size() != input_dim)
        {
          initial_param.u_list[i].setZero(input_dim);
        }
      }
    }
    Eigen::VectorXd planned_force_scales = ddp.planOnce(motion_param_func, ref_data_func, initial_param, t);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    const auto & motion_param = motion_param_func(t);
    const auto & ref_data = ref_data_func(t);
    const auto & planned_total_wrench = motion_param.calcTotalWrench(planned_force_scales, sim.state_.pos);
    ofs << t << " " << sim.state_.pos.transpose() << " " << ref_data.pos.transpose() << " "
        << sim.state_.vel.transpose() << " " << sim.state_.angular_momentum.transpose() << " "
        << planned_total_wrench.transpose() << " " << ddp.ddp_solver_->traceDataList().back().iter << " "
        << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((sim.state_.pos - ref_data.pos).norm(), 2.0); // [m]
    EXPECT_LT(sim.state_.vel.norm(), 2.0); // [m/s]
    EXPECT_LT(sim.state_.angular_momentum.norm(), 1.0); // [kg m^2/s]

    // Simulate
    t += sim_dt;
    sim.update(planned_total_wrench);

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
  EXPECT_LT((sim.state_.pos - ref_data.pos).norm(), 0.1); // [m]
  EXPECT_LT(sim.state_.vel.norm(), 0.1); // [m/s]
  EXPECT_LT(sim.state_.angular_momentum.norm(), 0.01); // [kg m^2/s]

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
            << "  plot \"" << file_path << "\" u 1:4 w lp, \"\" u 1:7 w l lw 2 # pos_z\n"
            << "  plot \"" << file_path << "\" u 1:12 w lp # angular_momentum_y\n"
            << "  plot \"" << file_path << "\" u 1:16 w lp # force_z\n"
            << "  plot \"" << file_path << "\" u 1:20 w lp # ddp_iter\n"
            << "  plot \"" << file_path << "\" u 1:21 w lp # computation_time\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
