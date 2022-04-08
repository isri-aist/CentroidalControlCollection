/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/DdpZmp.h>

using namespace CCC;

DdpZmp::DdpProblem::StateDimVector DdpZmp::DdpProblem::stateEq(double t,
                                                               const StateDimVector & x,
                                                               const InputDimVector & u) const
{
  StateDimVector x_dot;
  x_dot << x[1], (x[0] - u[0]) * u[2] / (mass_ * x[4]), x[3], (x[2] - u[1]) * u[2] / (mass_ * x[4]), x[5],
      u[2] / mass_ - constants::g;
  return x + dt_ * x_dot;
}

double DdpZmp::DdpProblem::runningCost(double t, const StateDimVector & x, const InputDimVector & u) const
{
  const auto & ref_data = ref_data_func_(t);

  return weight_param_.running_com_pos_z * 0.5 * std::pow(x[4] - ref_data.com_z, 2)
         + weight_param_.running_zmp * 0.5 * (u.head<2>() - ref_data.zmp).squaredNorm()
         + weight_param_.running_force_z * 0.5 * std::pow(u[2] - mass_ * constants::g, 2);
}

double DdpZmp::DdpProblem::terminalCost(double t, const StateDimVector & x) const
{
  const auto & ref_data = ref_data_func_(t);

  Eigen::Vector2d com_pos_xy;
  com_pos_xy << x[0], x[2];
  Eigen::Vector3d com_vel;
  com_vel << x[1], x[3], x[5];

  return weight_param_.terminal_com_pos_xy * 0.5 * (com_pos_xy - ref_data.zmp).squaredNorm()
         + weight_param_.terminal_com_pos_z * 0.5 * std::pow(x[4] - ref_data.com_z, 2)
         + weight_param_.terminal_com_vel * 0.5 * com_vel.squaredNorm();
}

void DdpZmp::DdpProblem::calcStateEqDeriv(double t,
                                          const StateDimVector & x,
                                          const InputDimVector & u,
                                          Eigen::Ref<StateStateDimMatrix> state_eq_deriv_x,
                                          Eigen::Ref<StateInputDimMatrix> state_eq_deriv_u) const
{
  state_eq_deriv_x.setZero();
  state_eq_deriv_x(0, 1) = 1;
  state_eq_deriv_x(1, 0) = u[2] / (mass_ * x[4]);
  state_eq_deriv_x(1, 4) = -1 * (x[0] - u[0]) * u[2] / (mass_ * std::pow(x[4], 2));
  state_eq_deriv_x(2, 3) = 1;
  state_eq_deriv_x(3, 2) = u[2] / (mass_ * x[4]);
  state_eq_deriv_x(3, 4) = -1 * (x[2] - u[1]) * u[2] / (mass_ * std::pow(x[4], 2));
  state_eq_deriv_x(4, 5) = 1;
  state_eq_deriv_x *= dt_;
  state_eq_deriv_x.diagonal().array() += 1.0;

  state_eq_deriv_u.setZero();
  state_eq_deriv_u(1, 0) = -1 * u[2] / (mass_ * x[4]);
  state_eq_deriv_u(1, 2) = (x[0] - u[0]) / (mass_ * x[4]);
  state_eq_deriv_u(3, 1) = -1 * u[2] / (mass_ * x[4]);
  state_eq_deriv_u(3, 2) = (x[2] - u[1]) / (mass_ * x[4]);
  state_eq_deriv_u(5, 2) = 1 / mass_;
  state_eq_deriv_u *= dt_;
}

void DdpZmp::DdpProblem::calcRunningCostDeriv(double t,
                                              const StateDimVector & x,
                                              const InputDimVector & u,
                                              Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                              Eigen::Ref<InputDimVector> running_cost_deriv_u) const
{
  const auto & ref_data = ref_data_func_(t);

  running_cost_deriv_x.setZero();
  running_cost_deriv_x[4] = weight_param_.running_com_pos_z * (x[4] - ref_data.com_z);

  running_cost_deriv_u << weight_param_.running_zmp * (u.head<2>() - ref_data.zmp),
      weight_param_.running_force_z * (u[2] - mass_ * constants::g);
}

void DdpZmp::DdpProblem::calcRunningCostDeriv(double t,
                                              const StateDimVector & x,
                                              const InputDimVector & u,
                                              Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                              Eigen::Ref<InputDimVector> running_cost_deriv_u,
                                              Eigen::Ref<StateStateDimMatrix> running_cost_deriv_xx,
                                              Eigen::Ref<InputInputDimMatrix> running_cost_deriv_uu,
                                              Eigen::Ref<StateInputDimMatrix> running_cost_deriv_xu) const
{
  const auto & ref_data = ref_data_func_(t);

  running_cost_deriv_x.setZero();
  running_cost_deriv_x[4] = weight_param_.running_com_pos_z * (x[4] - ref_data.com_z);

  running_cost_deriv_u << weight_param_.running_zmp * (u.head<2>() - ref_data.zmp),
      weight_param_.running_force_z * (u[2] - mass_ * constants::g);

  running_cost_deriv_xx.setZero();
  running_cost_deriv_xx.diagonal()[4] = weight_param_.running_com_pos_z;
  running_cost_deriv_uu.setZero();
  running_cost_deriv_uu.diagonal() << weight_param_.running_zmp, weight_param_.running_zmp,
      weight_param_.running_force_z;
  running_cost_deriv_xu.setZero();
}

void DdpZmp::DdpProblem::calcTerminalCostDeriv(double t,
                                               const StateDimVector & x,
                                               Eigen::Ref<StateDimVector> terminal_cost_deriv_x) const
{
  const auto & ref_data = ref_data_func_(t);

  terminal_cost_deriv_x[0] = weight_param_.terminal_com_pos_xy * (x[0] - ref_data.zmp[0]);
  terminal_cost_deriv_x[1] = weight_param_.terminal_com_vel * x[1];
  terminal_cost_deriv_x[2] = weight_param_.terminal_com_pos_xy * (x[2] - ref_data.zmp[1]);
  terminal_cost_deriv_x[3] = weight_param_.terminal_com_vel * x[3];
  terminal_cost_deriv_x[4] = weight_param_.terminal_com_pos_z * (x[4] - ref_data.com_z);
  terminal_cost_deriv_x[5] = weight_param_.terminal_com_vel * x[5];
}

void DdpZmp::DdpProblem::calcTerminalCostDeriv(double t,
                                               const StateDimVector & x,
                                               Eigen::Ref<StateDimVector> terminal_cost_deriv_x,
                                               Eigen::Ref<StateStateDimMatrix> terminal_cost_deriv_xx) const
{
  const auto & ref_data = ref_data_func_(t);

  terminal_cost_deriv_x[0] = weight_param_.terminal_com_pos_xy * (x[0] - ref_data.zmp[0]);
  terminal_cost_deriv_x[1] = weight_param_.terminal_com_vel * x[1];
  terminal_cost_deriv_x[2] = weight_param_.terminal_com_pos_xy * (x[2] - ref_data.zmp[1]);
  terminal_cost_deriv_x[3] = weight_param_.terminal_com_vel * x[3];
  terminal_cost_deriv_x[4] = weight_param_.terminal_com_pos_z * (x[4] - ref_data.com_z);
  terminal_cost_deriv_x[5] = weight_param_.terminal_com_vel * x[5];

  terminal_cost_deriv_xx.setZero();
  terminal_cost_deriv_xx.diagonal() << weight_param_.terminal_com_pos_xy, weight_param_.terminal_com_vel,
      weight_param_.terminal_com_pos_xy, weight_param_.terminal_com_vel, weight_param_.terminal_com_pos_z,
      weight_param_.terminal_com_vel;
}

DdpZmp::DdpProblem::StateDimVector DdpZmp::InitialParam::toState() const
{
  DdpProblem::StateDimVector x;
  x << pos[0], vel[0], pos[1], vel[1], pos[2], vel[2];
  return x;
}

DdpZmp::PlannedData DdpZmp::planOnce(const std::function<RefData(double)> & ref_data_func,
                                     const InitialParam & initial_param,
                                     double current_time,
                                     const std::vector<DdpProblem::InputDimVector> & initial_u_list)
{
  ddp_problem_->setRefDataFunc(ref_data_func);

  if(initial_u_list.empty())
  {
    ddp_solver_->solve(current_time, initial_param.toState(),
                       std::vector<DdpProblem::InputDimVector>(ddp_solver_->config().horizon_steps,
                                                               DdpProblem::InputDimVector::Zero()));
  }
  else
  {
    ddp_solver_->solve(current_time, initial_param.toState(), initial_u_list);
  }

  PlannedData planned_data;
  planned_data.zmp = ddp_solver_->controlData().u_list[0].head<2>();
  planned_data.force_z = ddp_solver_->controlData().u_list[0][2];

  return planned_data;
}
