/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DdpCentroidal.h>

#include "ContactManager.h"
#include "SimModels.h"

TEST(TestDdpCentroidal, PlanOnce)
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
  CCC::DdpCentroidal::WeightParam weight_param;
  weight_param.running_pos << 1.0, 1.0, 10.0;
  weight_param.terminal_pos << 1.0, 1.0, 10.0;
  CCC::DdpCentroidal ddp(mass, horizon_dt, horizon_steps, weight_param);
  ddp.ddp_solver_->config().initial_lambda = 1e-6;
  ddp.ddp_solver_->config().lambda_min = 1e-8;
  ddp.ddp_solver_->config().lambda_thre = 1e-7;

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
    ddp.ddp_solver_->config().max_iter = 1; // Set max_iter from second simulation iteration
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

TEST(TestDdpCentroidal, PlanLoop)
{
  double horizon_dt = 0.03; // [sec]
  double horizon_duration = 3.0; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.005; // [sec]
  double mass = 100.0; // [kg]

  // Setup DDP
  CCC::DdpCentroidal::WeightParam weight_param;
  weight_param.running_pos << 1.0, 1.0, 10.0;
  weight_param.terminal_pos << 1.0, 1.0, 10.0;
  CCC::DdpCentroidal ddp(mass, horizon_dt, horizon_steps, weight_param);
  ddp.ddp_solver_->config().initial_lambda = 1e-6;
  ddp.ddp_solver_->config().lambda_min = 1e-8;
  ddp.ddp_solver_->config().lambda_thre = 1e-7;

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
    return ref_data;
  };
  CCC::DdpCentroidal::InitialParam initial_param;
  initial_param.pos << 0.0, 0.0, 1.0; // [m]
  std::pair<double, double> motion_time_range(0.0, 3.0); // ([sec], [sec])

  // Plan
  ddp.planLoop(motion_param_func, ref_data_func, initial_param, motion_time_range, sim_dt);

  // Dump
  ddp.dumpMotionDataSeq("/tmp/TestDdpCentroidal_PlanLoop.txt", true);

  // Final check
  const auto & motion_data_final = ddp.motionDataSeq().rbegin()->second;
  const auto & ref_data = ref_data_func(motion_time_range.second);
  EXPECT_LT((motion_data_final.planned_pos - ref_data.pos).norm(), 0.1); // [m]
  EXPECT_LT(motion_data_final.planned_vel.norm(), 0.1); // [m/s]
  EXPECT_LT(motion_data_final.planned_angular_momentum.norm(), 0.01); // [kg m^2/s]
}

TEST(TestDdpCentroidal, CheckDerivatives)
{
  constexpr double deriv_eps = 1e-6;

  double horizon_dt = 0.03; // [sec]
  double mass = 100.0; // [Kg]
  CCC::DdpCentroidal::WeightParam weight_param;
  auto ddp_problem = std::make_shared<CCC::DdpCentroidal::DdpProblem>(horizon_dt, mass, weight_param);

  std::function<CCC::DdpCentroidal::MotionParam(double)> motion_param_func = [](double t) {
    CCC::DdpCentroidal::MotionParam motion_param;
    motion_param.vertex_ridge_list =
        makeVertexRidgeListFromRect({Eigen::Vector2d(-0.1, -0.1), Eigen::Vector2d(0.1, 0.1)});
    return motion_param;
  };
  std::function<CCC::DdpCentroidal::RefData(double)> ref_data_func = [](double t) {
    CCC::DdpCentroidal::RefData ref_data;
    ref_data.pos << 0.1, -0.2, 1.0; // [m]
    return ref_data;
  };
  ddp_problem->setMotionParamFunc(motion_param_func);
  ddp_problem->setRefDataFunc(ref_data_func);

  double t = 0;
  int state_dim = ddp_problem->stateDim();
  int input_dim = ddp_problem->inputDim(t);
  CCC::DdpCentroidal::DdpProblem::StateDimVector x;
  x << 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0;
  CCC::DdpCentroidal::DdpProblem::InputDimVector u(input_dim);
  u << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0;

  // Check state_eq_deriv
  {
    CCC::DdpCentroidal::DdpProblem::StateStateDimMatrix state_eq_deriv_x_analytical;
    CCC::DdpCentroidal::DdpProblem::StateInputDimMatrix state_eq_deriv_u_analytical(state_dim, input_dim);
    ddp_problem->calcStateEqDeriv(t, x, u, state_eq_deriv_x_analytical, state_eq_deriv_u_analytical);

    CCC::DdpCentroidal::DdpProblem::StateStateDimMatrix state_eq_deriv_x_numerical;
    CCC::DdpCentroidal::DdpProblem::StateInputDimMatrix state_eq_deriv_u_numerical(state_dim, input_dim);
    for(int i = 0; i < state_dim; i++)
    {
      state_eq_deriv_x_numerical.col(i) =
          (ddp_problem->stateEq(t, x + deriv_eps * CCC::DdpCentroidal::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->stateEq(t, x - deriv_eps * CCC::DdpCentroidal::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < input_dim; i++)
    {
      state_eq_deriv_u_numerical.col(i) =
          (ddp_problem->stateEq(t, x,
                                u + deriv_eps * CCC::DdpCentroidal::DdpProblem::InputDimVector::Unit(input_dim, i))
           - ddp_problem->stateEq(t, x,
                                  u - deriv_eps * CCC::DdpCentroidal::DdpProblem::InputDimVector::Unit(input_dim, i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((state_eq_deriv_x_analytical - state_eq_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((state_eq_deriv_u_analytical - state_eq_deriv_u_numerical).norm(), 1e-6);
  }

  // Check running_cost_deriv
  {
    CCC::DdpCentroidal::DdpProblem::StateDimVector running_cost_deriv_x_analytical;
    CCC::DdpCentroidal::DdpProblem::InputDimVector running_cost_deriv_u_analytical(input_dim);
    ddp_problem->calcRunningCostDeriv(t, x, u, running_cost_deriv_x_analytical, running_cost_deriv_u_analytical);

    CCC::DdpCentroidal::DdpProblem::StateDimVector running_cost_deriv_x_numerical;
    CCC::DdpCentroidal::DdpProblem::InputDimVector running_cost_deriv_u_numerical(input_dim);
    for(int i = 0; i < state_dim; i++)
    {
      running_cost_deriv_x_numerical[i] =
          (ddp_problem->runningCost(t, x + deriv_eps * CCC::DdpCentroidal::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->runningCost(t, x - deriv_eps * CCC::DdpCentroidal::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < input_dim; i++)
    {
      running_cost_deriv_u_numerical[i] =
          (ddp_problem->runningCost(t, x,
                                    u + deriv_eps * CCC::DdpCentroidal::DdpProblem::InputDimVector::Unit(input_dim, i))
           - ddp_problem->runningCost(
               t, x, u - deriv_eps * CCC::DdpCentroidal::DdpProblem::InputDimVector::Unit(input_dim, i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((running_cost_deriv_x_analytical - running_cost_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((running_cost_deriv_u_analytical - running_cost_deriv_u_numerical).norm(), 1e-6);
  }

  // Check terminal_cost_deriv
  {
    CCC::DdpCentroidal::DdpProblem::StateDimVector terminal_cost_deriv_x_analytical;
    ddp_problem->calcTerminalCostDeriv(t, x, terminal_cost_deriv_x_analytical);

    CCC::DdpCentroidal::DdpProblem::StateDimVector terminal_cost_deriv_x_numerical;
    for(int i = 0; i < state_dim; i++)
    {
      terminal_cost_deriv_x_numerical[i] =
          (ddp_problem->terminalCost(t, x + deriv_eps * CCC::DdpCentroidal::DdpProblem::StateDimVector::Unit(i))
           - ddp_problem->terminalCost(t, x - deriv_eps * CCC::DdpCentroidal::DdpProblem::StateDimVector::Unit(i)))
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
