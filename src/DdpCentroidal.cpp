/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/DdpCentroidal.h>

using namespace CCC;

namespace
{
/** \brief Calculate a matrix corresponding to the cross product. */
Eigen::Matrix3d crossMat(const Eigen::Vector3d & vec)
{
  Eigen::Matrix3d mat;
  mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
  return mat;
}
} // namespace

Eigen::Vector6d DdpCentroidal::MotionParam::calcTotalWrench(const Eigen::VectorXd & force_scales,
                                                            const Eigen::Vector3d & moment_origin) const
{
  const Eigen::Ref<const Eigen::Matrix3Xd> & vertices_mat = vertex_ridge_list.topRows<3>();
  const Eigen::Ref<const Eigen::Matrix3Xd> & ridges_mat = vertex_ridge_list.bottomRows<3>();

  Eigen::Vector6d wrench;
  wrench.head<3>() = ridges_mat * force_scales;
  wrench.tail<3>().setZero();
  for(int i = 0; i < force_scales.size(); i++)
  {
    wrench.tail<3>() += force_scales[i] * (vertices_mat.col(i) - moment_origin).cross(ridges_mat.col(i));
  }
  return wrench;
}

DdpCentroidal::DdpProblem::StateDimVector DdpCentroidal::DdpProblem::stateEq(double t,
                                                                             const StateDimVector & x,
                                                                             const InputDimVector & u) const
{
  const MotionParam & motion_param = motion_param_func_(t);
  const Eigen::Ref<const Eigen::Matrix3Xd> & vertices_mat = motion_param.vertex_ridge_list.topRows<3>();
  const Eigen::Ref<const Eigen::Matrix3Xd> & ridges_mat = motion_param.vertex_ridge_list.bottomRows<3>();

  const Eigen::Ref<const Eigen::Vector3d> & pos = x.segment<3>(0);
  const Eigen::Ref<const Eigen::Vector3d> & linear_momentum = x.segment<3>(3);
  const Eigen::Ref<const Eigen::Vector3d> & angular_momentum = x.segment<3>(6);

  StateDimVector x_dot;
  Eigen::Ref<Eigen::Vector3d> pos_dot = x_dot.segment<3>(0);
  Eigen::Ref<Eigen::Vector3d> linear_momentum_dot = x_dot.segment<3>(3);
  Eigen::Ref<Eigen::Vector3d> angular_momentum_dot = x_dot.segment<3>(6);
  pos_dot = linear_momentum / mass_;
  linear_momentum_dot = ridges_mat * u - mass_ * Eigen::Vector3d(0, 0, constants::g);
  angular_momentum_dot.setZero();
  for(int i = 0; i < u.size(); i++)
  {
    angular_momentum_dot += u[i] * (vertices_mat.col(i) - pos).cross(ridges_mat.col(i));
  }

  return x + dt_ * x_dot;
}

double DdpCentroidal::DdpProblem::runningCost(double t, const StateDimVector & x, const InputDimVector & u) const
{
  const auto & ref_data = ref_data_func_(t);

  return 0.5 * weight_param_.running_pos.dot((x.segment<3>(0) - ref_data.pos).cwiseAbs2())
         + 0.5 * weight_param_.running_linear_momentum.dot((x.segment<3>(3)).cwiseAbs2())
         + 0.5 * weight_param_.running_angular_momentum.dot((x.segment<3>(6)).cwiseAbs2())
         + 0.5 * weight_param_.running_force * u.squaredNorm();
}

double DdpCentroidal::DdpProblem::terminalCost(double t, const StateDimVector & x) const
{
  const auto & ref_data = ref_data_func_(t);

  return 0.5 * weight_param_.terminal_pos.dot((x.segment<3>(0) - ref_data.pos).cwiseAbs2())
         + 0.5 * weight_param_.terminal_linear_momentum.dot((x.segment<3>(3)).cwiseAbs2())
         + 0.5 * weight_param_.terminal_angular_momentum.dot((x.segment<3>(6)).cwiseAbs2());
}

void DdpCentroidal::DdpProblem::calcStateEqDeriv(double t,
                                                 const StateDimVector & x,
                                                 const InputDimVector & u,
                                                 Eigen::Ref<StateStateDimMatrix> state_eq_deriv_x,
                                                 Eigen::Ref<StateInputDimMatrix> state_eq_deriv_u) const
{
  const MotionParam & motion_param = motion_param_func_(t);
  const Eigen::Ref<const Eigen::Matrix3Xd> & vertices_mat = motion_param.vertex_ridge_list.topRows<3>();
  const Eigen::Ref<const Eigen::Matrix3Xd> & ridges_mat = motion_param.vertex_ridge_list.bottomRows<3>();

  const Eigen::Ref<const Eigen::Vector3d> & pos = x.segment<3>(0);
  const Eigen::Ref<const Eigen::Vector3d> & linear_momentum = x.segment<3>(3);
  const Eigen::Ref<const Eigen::Vector3d> & angular_momentum = x.segment<3>(6);

  state_eq_deriv_x.setZero();
  state_eq_deriv_x.block<3, 3>(0, 3).diagonal().setConstant(1 / mass_);
  state_eq_deriv_x.block<3, 3>(6, 0) = crossMat(ridges_mat * u);
  state_eq_deriv_x *= dt_;
  state_eq_deriv_x += StateStateDimMatrix::Identity();

  state_eq_deriv_u.setZero();
  state_eq_deriv_u.middleRows<3>(3) = ridges_mat;
  for(int i = 0; i < u.size(); i++)
  {
    state_eq_deriv_u.middleRows<3>(6).col(i) = (vertices_mat.col(i) - pos).cross(ridges_mat.col(i));
  }
  state_eq_deriv_u *= dt_;
}

void DdpCentroidal::DdpProblem::calcRunningCostDeriv(double t,
                                                     const StateDimVector & x,
                                                     const InputDimVector & u,
                                                     Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                                     Eigen::Ref<InputDimVector> running_cost_deriv_u) const
{
  const auto & ref_data = ref_data_func_(t);

  running_cost_deriv_x << weight_param_.running_pos.cwiseProduct(x.segment<3>(0) - ref_data.pos),
      weight_param_.running_linear_momentum.cwiseProduct(x.segment<3>(3)),
      weight_param_.running_angular_momentum.cwiseProduct(x.segment<3>(6));
  running_cost_deriv_u = weight_param_.running_force * u;
}

void DdpCentroidal::DdpProblem::calcRunningCostDeriv(double t,
                                                     const StateDimVector & x,
                                                     const InputDimVector & u,
                                                     Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                                     Eigen::Ref<InputDimVector> running_cost_deriv_u,
                                                     Eigen::Ref<StateStateDimMatrix> running_cost_deriv_xx,
                                                     Eigen::Ref<InputInputDimMatrix> running_cost_deriv_uu,
                                                     Eigen::Ref<StateInputDimMatrix> running_cost_deriv_xu) const
{
  calcRunningCostDeriv(t, x, u, running_cost_deriv_x, running_cost_deriv_u);

  running_cost_deriv_xx.setZero();
  running_cost_deriv_xx.diagonal() << weight_param_.running_pos, weight_param_.running_linear_momentum,
      weight_param_.running_angular_momentum;
  running_cost_deriv_uu.setIdentity();
  running_cost_deriv_uu *= weight_param_.running_force;
  running_cost_deriv_xu.setZero();
}

void DdpCentroidal::DdpProblem::calcTerminalCostDeriv(double t,
                                                      const StateDimVector & x,
                                                      Eigen::Ref<StateDimVector> terminal_cost_deriv_x) const
{
  const auto & ref_data = ref_data_func_(t);

  terminal_cost_deriv_x << weight_param_.terminal_pos.cwiseProduct(x.segment<3>(0) - ref_data.pos),
      weight_param_.terminal_linear_momentum.cwiseProduct(x.segment<3>(3)),
      weight_param_.terminal_angular_momentum.cwiseProduct(x.segment<3>(6));
}

void DdpCentroidal::DdpProblem::calcTerminalCostDeriv(double t,
                                                      const StateDimVector & x,
                                                      Eigen::Ref<StateDimVector> terminal_cost_deriv_x,
                                                      Eigen::Ref<StateStateDimMatrix> terminal_cost_deriv_xx) const
{
  calcTerminalCostDeriv(t, x, terminal_cost_deriv_x);

  terminal_cost_deriv_xx.setZero();
  terminal_cost_deriv_xx.diagonal() << weight_param_.terminal_pos, weight_param_.terminal_linear_momentum,
      weight_param_.terminal_angular_momentum;
}

DdpCentroidal::DdpProblem::StateDimVector DdpCentroidal::InitialParam::toState(double mass) const
{
  DdpProblem::StateDimVector state;
  state << pos, mass * vel, angular_momentum;
  return state;
}

DdpCentroidal::DdpCentroidal(double mass, double horizon_dt, int horizon_steps, const WeightParam & weight_param)
: ddp_problem_(std::make_shared<DdpProblem>(horizon_dt, mass, weight_param)),
  ddp_solver_(std::make_shared<nmpc_ddp::DDPSolver<9, Eigen::Dynamic>>(ddp_problem_))
{
  ddp_solver_->config().horizon_steps = horizon_steps;
}

Eigen::VectorXd DdpCentroidal::planOnce(const std::function<MotionParam(double)> & motion_param_func,
                                        const std::function<RefData(double)> & ref_data_func,
                                        const InitialParam & initial_param,
                                        double current_time)
{
  ddp_problem_->setMotionParamFunc(motion_param_func);
  ddp_problem_->setRefDataFunc(ref_data_func);

  if(initial_param.u_list.empty())
  {
    std::vector<DdpProblem::InputDimVector> current_u_list;
    for(int i = 0; i < ddp_solver_->config().horizon_steps; i++)
    {
      double t = current_time + i * ddp_problem_->dt();
      current_u_list.push_back(DdpProblem::InputDimVector::Zero(ddp_problem_->inputDim(t)));
    }
    ddp_solver_->solve(current_time, initial_param.toState(ddp_problem_->mass_), current_u_list);
  }
  else
  {
    ddp_solver_->solve(current_time, initial_param.toState(ddp_problem_->mass_), initial_param.u_list);
  }

  return ddp_solver_->controlData().u_list[0];
}
