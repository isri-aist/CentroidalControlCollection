/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DdpZmp.h>
#include <CCC/EigenTypes.h>

#include "SimModels.h"

/** \brief CoM position and velocity. */
struct ComPosVel
{
  Eigen::Vector2d x = Eigen::Vector2d::Zero();
  Eigen::Vector2d y = Eigen::Vector2d::Zero();
  Eigen::Vector2d z = Eigen::Vector2d::Zero();
};

TEST(TestDdpZmp, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.02; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.005; // [sec]
  double mass = 100.0; // [Kg]
  double ref_com_height = 1.0; // [m]

  // Setup DDP
  CCC::DdpZmp ddp(mass, horizon_dt, horizon_steps);
  ddp.ddp_solver_->config().max_iter = 3;

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
  std::function<CCC::DdpZmp::RefData(double)> ref_data_func = [&](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::DdpZmp::RefData ref_data;
    ref_data.zmp.setConstant(ref_zmp_func(t));
    ref_data.com_z = ref_com_height;

    return ref_data;
  };

  // Setup dump file
  std::string file_path = "/tmp/TestDdpZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos com_vel planned_zmp ref_zmp ddp_iter" << std::endl;

  // Setup control loop
  ComPosVel com_pos_vel;
  com_pos_vel.z[0] = ref_com_height;
  CCC::DdpZmp::PlannedData planned_data;

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  bool first_iter = true;
  while(t < end_time)
  {
    // Plan
    CCC::DdpZmp::InitialParam initial_param;
    initial_param.pos << com_pos_vel.x[0], com_pos_vel.y[0], com_pos_vel.z[0];
    initial_param.vel << com_pos_vel.x[1], com_pos_vel.y[1], com_pos_vel.z[1];
    if(first_iter)
    {
      first_iter = false;
      initial_param.u_list.assign(horizon_steps, CCC::DdpZmp::DdpProblem::InputDimVector(
                                                     com_pos_vel.x[0], com_pos_vel.y[0], mass * CCC::constants::g));
    }
    else
    {
      initial_param.u_list = ddp.ddp_solver_->controlData().u_list;
    }
    planned_data = ddp.planOnce(ref_data_func, initial_param, t);

    // Dump
    const auto & ref_data = ref_data_func(t);
    ofs << t << " " << com_pos_vel.x.transpose() << " " << planned_data.zmp.x() << " " << ref_data.zmp.x() << " "
        << ddp.ddp_solver_->traceDataList().back().iter << std::endl;

    // Check
    EXPECT_LT(std::abs(planned_data.zmp.x() - ref_data.zmp.x()), 0.1); // [m]

    // Setup simulation
    ComZmpSimModel1d sim_model_xy(com_pos_vel.z[0]);
    sim_model_xy.calcDiscMatrix(sim_dt);
    VerticalSimModel sim_model_z(mass);
    sim_model_z.calcDiscMatrix(sim_dt);

    // Simulate
    t += sim_dt;
    Eigen::Vector1d sim_input;
    sim_input << planned_data.zmp.x();
    com_pos_vel.x = sim_model_xy.stateEqDisc(com_pos_vel.x, sim_input);
    sim_input << planned_data.zmp.y();
    com_pos_vel.y = sim_model_xy.stateEqDisc(com_pos_vel.y, sim_input);
    sim_input << planned_data.force_z;
    com_pos_vel.z = sim_model_z.stateEqDisc(com_pos_vel.z, sim_input);
  }

  // Final check
  const auto & ref_data = ref_data_func(t);
  EXPECT_LT(std::abs(planned_data.zmp.x() - ref_data.zmp.x()), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel.x[0] - ref_data.zmp.x()), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel.x[1]), 1e-2);

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:5 w l lw 2\n";
}

TEST(TestDdpZmp, CheckDerivatives)
{
  constexpr double deriv_eps = 1e-6;

  double horizon_dt = 0.005; // [sec]
  double mass = 100.0; // [Kg]
  CCC::DdpZmp::WeightParam weight_param;
  auto ddp_problem = std::make_shared<CCC::DdpZmp::DdpProblem>(horizon_dt, mass, weight_param);

  std::function<CCC::DdpZmp::RefData(double)> ref_data_func = [](double t) {
    CCC::DdpZmp::RefData ref_data;
    ref_data.zmp << 0.1, -0.2; // [m]
    ref_data.com_z = 1.0; // [m]
    return ref_data;
  };
  ddp_problem->setRefDataFunc(ref_data_func);

  double t = 0;
  CCC::DdpZmp::DdpProblem::StateDimVector x;
  x << 1.0, 0.1, -2.1, -0.5, 1.1, 0.5; // [m], [m/s]
  CCC::DdpZmp::DdpProblem::InputDimVector u;
  u << 1.0, -2.0, 1000.0; // [m], [N]

  // Check state_eq_deriv
  {
    CCC::DdpZmp::DdpProblem::StateStateDimMatrix state_eq_deriv_x_analytical;
    CCC::DdpZmp::DdpProblem::StateInputDimMatrix state_eq_deriv_u_analytical;
    ddp_problem->calcStateEqDeriv(t, x, u, state_eq_deriv_x_analytical, state_eq_deriv_u_analytical);

    CCC::DdpZmp::DdpProblem::StateStateDimMatrix state_eq_deriv_x_numerical;
    CCC::DdpZmp::DdpProblem::StateInputDimMatrix state_eq_deriv_u_numerical;
    for(int i = 0; i < ddp_problem->stateDim(); i++)
    {
      state_eq_deriv_x_numerical.col(i) =
          (ddp_problem->stateEq(t, x + deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->stateEq(t, x - deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < ddp_problem->inputDim(); i++)
    {
      state_eq_deriv_u_numerical.col(i) =
          (ddp_problem->stateEq(t, x, u + deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i))
           - ddp_problem->stateEq(t, x, u - deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((state_eq_deriv_x_analytical - state_eq_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((state_eq_deriv_u_analytical - state_eq_deriv_u_numerical).norm(), 1e-6);
  }

  // Check running_cost_deriv
  {
    CCC::DdpZmp::DdpProblem::StateDimVector running_cost_deriv_x_analytical;
    CCC::DdpZmp::DdpProblem::InputDimVector running_cost_deriv_u_analytical;
    ddp_problem->calcRunningCostDeriv(t, x, u, running_cost_deriv_x_analytical, running_cost_deriv_u_analytical);

    CCC::DdpZmp::DdpProblem::StateDimVector running_cost_deriv_x_numerical;
    CCC::DdpZmp::DdpProblem::InputDimVector running_cost_deriv_u_numerical;
    for(int i = 0; i < ddp_problem->stateDim(); i++)
    {
      running_cost_deriv_x_numerical[i] =
          (ddp_problem->runningCost(t, x + deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->runningCost(t, x - deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < ddp_problem->inputDim(); i++)
    {
      running_cost_deriv_u_numerical[i] =
          (ddp_problem->runningCost(t, x, u + deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i))
           - ddp_problem->runningCost(t, x, u - deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((running_cost_deriv_x_analytical - running_cost_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((running_cost_deriv_u_analytical - running_cost_deriv_u_numerical).norm(), 1e-6);
  }

  // Check terminal_cost_deriv
  {
    CCC::DdpZmp::DdpProblem::StateDimVector terminal_cost_deriv_x_analytical;
    ddp_problem->calcTerminalCostDeriv(t, x, terminal_cost_deriv_x_analytical);

    CCC::DdpZmp::DdpProblem::StateDimVector terminal_cost_deriv_x_numerical;
    for(int i = 0; i < ddp_problem->stateDim(); i++)
    {
      terminal_cost_deriv_x_numerical[i] =
          (ddp_problem->terminalCost(t, x + deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i))
           - ddp_problem->terminalCost(t, x - deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i)))
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
