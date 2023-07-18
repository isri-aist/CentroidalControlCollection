/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DdpSingleRigidBody.h>

#include "ContactManager.h"
#include "SimModels.h"

TEST(TestDdpSingleRigidBody, PlanOnce)
{
  double horizon_dt = 0.03; // [sec]
  double horizon_duration = 3.0; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.005; // [sec]
  double mass = 100.0; // [kg]
  Eigen::Vector3d moment_of_inertia = Eigen::Vector3d(40.0, 20.0, 10.0); // [kg m^2]
  std::vector<double> disturb_time_list = {1.0}; // [sec]
  sva::ForceVecd disturb_impulse_per_mass =
      sva::ForceVecd(Eigen::Vector3d::Zero(), Eigen::Vector3d(0.05, 0.05, 0.0)); // [m/s], [m^2/s]

  // Setup DDP
  std::vector<double> computation_duration_list;
  CCC::DdpSingleRigidBody::WeightParam weight_param;
  weight_param.running_pos << 1.0, 1.0, 10.0;
  weight_param.running_ori << 0.5, 0.5, 0.5;
  weight_param.terminal_pos << 1.0, 1.0, 10.0;
  weight_param.terminal_ori << 0.5, 0.5, 0.5;
  CCC::DdpSingleRigidBody ddp(mass, horizon_dt, horizon_steps, weight_param);

  // Setup contact
  std::function<CCC::DdpSingleRigidBody::MotionParam(double)> motion_param_func = [moment_of_inertia](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::DdpSingleRigidBody::MotionParam motion_param;
    if(t < 1.4)
    {
      motion_param.contact_list.push_back(
          makeContactFromRect({Eigen::Vector2d(-0.1, -0.5), Eigen::Vector2d(0.1, 0.5)}));
    }
    else if(t < 1.6)
    {
    }
    else
    {
      motion_param.contact_list.push_back(makeContactFromRect({Eigen::Vector2d(0.4, -0.5), Eigen::Vector2d(0.6, 0.5)}));
    }
    motion_param.inertia_mat.diagonal() = moment_of_inertia;
    return motion_param;
  };
  std::function<CCC::DdpSingleRigidBody::RefData(double)> ref_data_func = [](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::DdpSingleRigidBody::RefData ref_data;
    if(t < 1.4)
    {
      ref_data.pos << 0.0, 0.0, 1.0; // [m]
    }
    else if(t < 1.6)
    {
      ref_data.pos << 0.25, 0.0, 1.2; // [m]
    }
    else
    {
      ref_data.pos << 0.5, 0.0, 1.0; // [m]
    }
    if(2.2 < t && t < 2.4)
    {
      ref_data.ori << 0.0, 0.0, 0.3; // [rad]
    }
    else
    {
      ref_data.ori << 0.0, 0.0, 0.0; // [rad]
    }
    return ref_data;
  };

  // Setup simulation
  CentroidalSim sim(mass, moment_of_inertia, sim_dt);
  sim.state_.pos.linear() = ref_data_func(0.0).pos;
  sim.state_.pos.angular() = ref_data_func(0.0).ori.reverse();

  // Setup dump file
  std::string file_path = "/tmp/TestDdpSingleRigidBody.txt";
  std::ofstream ofs(file_path);
  ofs << "time planned_com_pos_x planned_com_pos_y planned_com_pos_z ref_com_pos_x ref_com_pos_y ref_com_pos_z "
         "planned_base_ori_x planned_base_ori_y planned_base_ori_z ref_base_ori_x ref_base_ori_y ref_base_ori_z "
         "planned_linear_vel_x planned_linear_vel_y planned_linear_vel_z "
         "planned_angular_vel_x planned_angular_vel_y planned_angular_vel_z "
         "wrench_nx wrench_ny wrench_nz wrench_fx wrench_fy wrench_fz ddp_iter "
         "computation_time"
      << std::endl;

  // Run control loop
  constexpr double end_time = 3.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    CCC::DdpSingleRigidBody::InitialParam initial_param;
    initial_param.pos = sim.state_.pos.linear();
    // Note that sim.state_.pos is in the order (X,Y,Z), but initial_param.ori is in the order (Z,Y,X).
    initial_param.ori = sim.state_.pos.angular().reverse();
    initial_param.linear_vel = sim.state_.vel.linear();
    initial_param.angular_vel = sim.state_.vel.angular();
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
    ddp.ddp_solver_->config().max_iter = 1; // Set max_iter from second simulation iteration
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    const auto & motion_param = motion_param_func(t);
    const auto & ref_data = ref_data_func(t);
    const auto & planned_total_wrench =
        ForceColl::calcTotalWrench(motion_param.contact_list, planned_force_scales, sim.state_.pos.linear());
    ofs << t << " " << sim.state_.pos.linear().transpose() << " " << ref_data.pos.transpose() << " "
        << sim.state_.pos.angular().transpose() << " " << ref_data.ori.reverse().transpose() << " "
        << sim.state_.vel.linear().transpose() << " " << sim.state_.vel.angular().transpose() << " "
        << planned_total_wrench.vector().transpose() << " " << ddp.ddp_solver_->traceDataList().back().iter << " "
        << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((sim.state_.pos.linear() - ref_data.pos).norm(), 2.0); // [m]
    EXPECT_LT((sim.state_.pos.angular() - ref_data.ori).norm(), 1.0); // [rad]
    EXPECT_LT(sim.state_.vel.linear().norm(), 2.0); // [m/s]
    EXPECT_LT(sim.state_.vel.angular().norm(), 2.0); // [rad/s]

    // Simulate
    t += sim_dt;
    sim.update(CentroidalSim::Input(planned_total_wrench));

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
  EXPECT_LT((sim.state_.pos.linear() - ref_data.pos).norm(), 0.1); // [m]
  EXPECT_LT((sim.state_.pos.angular() - ref_data.ori).norm(), 0.1); // [rad]
  EXPECT_LT(sim.state_.vel.linear().norm(), 0.1); // [m/s]
  EXPECT_LT(sim.state_.vel.angular().norm(), 0.1); // [rad/s]

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
            << "  plot \"" << file_path << "\" u 1:8 w lp, \"\" u 1:11 w l lw 2 # ori_x\n"
            << "  plot \"" << file_path << "\" u 1:9 w lp, \"\" u 1:12 w l lw 2 # ori_y\n"
            << "  plot \"" << file_path << "\" u 1:20 w lp # moment_x\n"
            << "  plot \"" << file_path << "\" u 1:25 w lp # force_z\n"
            << "  plot \"" << file_path << "\" u 1:26 w lp # ddp_iter\n"
            << "  plot \"" << file_path << "\" u 1:27 w lp # computation_time\n";
}

TEST(TestDdpSingleRigidBody, CheckDerivatives)
{
  constexpr double deriv_eps = 1e-6;

  double horizon_dt = 0.03; // [sec]
  double mass = 100.0; // [Kg]
  CCC::DdpSingleRigidBody::WeightParam weight_param;
  auto ddp_problem = std::make_shared<CCC::DdpSingleRigidBody::DdpProblem>(horizon_dt, mass, weight_param);

  std::function<CCC::DdpSingleRigidBody::MotionParam(double)> motion_param_func = [](double // t
                                                                                  ) {
    CCC::DdpSingleRigidBody::MotionParam motion_param;
    motion_param.contact_list.push_back(makeContactFromRect({Eigen::Vector2d(-0.1, -0.1), Eigen::Vector2d(0.1, 0.1)}));
    motion_param.inertia_mat.diagonal() << 15.0, 10.0, 5.0;
    return motion_param;
  };
  std::function<CCC::DdpSingleRigidBody::RefData(double)> ref_data_func = [](double // t
                                                                          ) {
    CCC::DdpSingleRigidBody::RefData ref_data;
    ref_data.pos << 0.1, -0.2, 1.0; // [m]
    ref_data.ori << -0.1, 0.2, -0.3; // [rad]
    return ref_data;
  };
  ddp_problem->setMotionParamFunc(motion_param_func);
  ddp_problem->setRefDataFunc(ref_data_func);

  double t = 0;
  int state_dim = ddp_problem->stateDim();
  int input_dim = ddp_problem->inputDim(t);
  CCC::DdpSingleRigidBody::DdpProblem::StateDimVector x;
  x << 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0, -10.0, 11.0, -12.0;
  CCC::DdpSingleRigidBody::DdpProblem::InputDimVector u(input_dim);
  u << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0;

  // Check state_eq_deriv
  {
    CCC::DdpSingleRigidBody::DdpProblem::StateStateDimMatrix state_eq_deriv_x_analytical;
    CCC::DdpSingleRigidBody::DdpProblem::StateInputDimMatrix state_eq_deriv_u_analytical(state_dim, input_dim);
    ddp_problem->calcStateEqDeriv(t, x, u, state_eq_deriv_x_analytical, state_eq_deriv_u_analytical);

    CCC::DdpSingleRigidBody::DdpProblem::StateStateDimMatrix state_eq_deriv_x_numerical;
    CCC::DdpSingleRigidBody::DdpProblem::StateInputDimMatrix state_eq_deriv_u_numerical(state_dim, input_dim);
    for(int i = 0; i < state_dim; i++)
    {
      state_eq_deriv_x_numerical.col(i) =
          (ddp_problem->stateEq(t, x + deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->stateEq(t, x - deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < input_dim; i++)
    {
      state_eq_deriv_u_numerical.col(i) =
          (ddp_problem->stateEq(t, x,
                                u + deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::InputDimVector::Unit(input_dim, i))
           - ddp_problem->stateEq(
               t, x, u - deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::InputDimVector::Unit(input_dim, i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((state_eq_deriv_x_analytical - state_eq_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((state_eq_deriv_u_analytical - state_eq_deriv_u_numerical).norm(), 1e-6);
  }

  // Check running_cost_deriv
  {
    CCC::DdpSingleRigidBody::DdpProblem::StateDimVector running_cost_deriv_x_analytical;
    CCC::DdpSingleRigidBody::DdpProblem::InputDimVector running_cost_deriv_u_analytical(input_dim);
    ddp_problem->calcRunningCostDeriv(t, x, u, running_cost_deriv_x_analytical, running_cost_deriv_u_analytical);

    CCC::DdpSingleRigidBody::DdpProblem::StateDimVector running_cost_deriv_x_numerical;
    CCC::DdpSingleRigidBody::DdpProblem::InputDimVector running_cost_deriv_u_numerical(input_dim);
    for(int i = 0; i < state_dim; i++)
    {
      running_cost_deriv_x_numerical[i] =
          (ddp_problem->runningCost(t, x + deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->runningCost(t, x - deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::StateDimVector::Unit(i),
                                      u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < input_dim; i++)
    {
      running_cost_deriv_u_numerical[i] =
          (ddp_problem->runningCost(
               t, x, u + deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::InputDimVector::Unit(input_dim, i))
           - ddp_problem->runningCost(
               t, x, u - deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::InputDimVector::Unit(input_dim, i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((running_cost_deriv_x_analytical - running_cost_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((running_cost_deriv_u_analytical - running_cost_deriv_u_numerical).norm(), 1e-6);
  }

  // Check terminal_cost_deriv
  {
    CCC::DdpSingleRigidBody::DdpProblem::StateDimVector terminal_cost_deriv_x_analytical;
    ddp_problem->calcTerminalCostDeriv(t, x, terminal_cost_deriv_x_analytical);

    CCC::DdpSingleRigidBody::DdpProblem::StateDimVector terminal_cost_deriv_x_numerical;
    for(int i = 0; i < state_dim; i++)
    {
      terminal_cost_deriv_x_numerical[i] =
          (ddp_problem->terminalCost(t, x + deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::StateDimVector::Unit(i))
           - ddp_problem->terminalCost(t, x - deriv_eps * CCC::DdpSingleRigidBody::DdpProblem::StateDimVector::Unit(i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((terminal_cost_deriv_x_analytical - terminal_cost_deriv_x_numerical).norm(), 1e-6);
  }
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
