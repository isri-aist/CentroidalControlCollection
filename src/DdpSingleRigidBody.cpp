/* Author: Masaki Murooka */

#include <ForceColl/Contact.h>

#include <CCC/Constants.h>
#include <CCC/DdpSingleRigidBody.h>

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

/** \brief Calculate a matrix to convert the angular velocity to the time derivative of the Euler angles.
    \param ori Euler angles

    ZYX Euler angles are assumed. The 3D vector representing the Euler angles is assumed to contain the rotation angles
   around the Z, Y, and X axes in this order.
 */
Eigen::Matrix3d matAngularVelToEulerDot(const Eigen::Vector3d & ori)
{
  // See "Time Derivative of Euler Angles ZYX <=> Angular Velocity" in
  // https://opensource.docs.anymal.com/doxygen/kindr/master/cheatsheet_latest.pdf
  double cos_alpha = std::cos(ori.x());
  double sin_alpha = std::sin(ori.x());
  double cos_beta = std::cos(ori.y());
  double sin_beta = std::sin(ori.y());
  Eigen::Matrix3d euler_trans_mat;
  euler_trans_mat << (cos_alpha * sin_beta) / cos_beta, (sin_beta * sin_alpha) / cos_beta, 1.0, -1 * sin_alpha,
      cos_alpha, 0.0, cos_alpha / cos_beta, sin_alpha / cos_beta, 0.0;
  return euler_trans_mat;
}
} // namespace

int DdpSingleRigidBody::DdpProblem::inputDim(double t) const
{
  const MotionParam & motion_param = motion_param_func_(t);
  int input_dim = 0;
  for(const auto & contact : motion_param.contact_list)
  {
    input_dim += contact->ridgeNum();
  }
  return input_dim;
}

DdpSingleRigidBody::DdpProblem::StateDimVector DdpSingleRigidBody::DdpProblem::stateEq(double t,
                                                                                       const StateDimVector & x,
                                                                                       const InputDimVector & u) const
{
  const MotionParam & motion_param = motion_param_func_(t);
  const Eigen::Ref<const Eigen::Matrix3d> inertia_mat = motion_param.inertia_mat;

  const Eigen::Ref<const Eigen::Vector3d> & pos = x.segment<3>(0);
  const Eigen::Ref<const Eigen::Vector3d> & ori = x.segment<3>(3);
  const Eigen::Ref<const Eigen::Vector3d> & linear_vel = x.segment<3>(6);
  const Eigen::Ref<const Eigen::Vector3d> & angular_vel = x.segment<3>(9);

  StateDimVector x_dot;
  Eigen::Ref<Eigen::Vector3d> pos_dot = x_dot.segment<3>(0);
  Eigen::Ref<Eigen::Vector3d> ori_dot = x_dot.segment<3>(3);
  Eigen::Ref<Eigen::Vector3d> linear_vel_dot = x_dot.segment<3>(6);
  Eigen::Ref<Eigen::Vector3d> angular_vel_dot = x_dot.segment<3>(9);
  pos_dot = linear_vel;
  // \todo Support other Euler angles
  ori_dot = matAngularVelToEulerDot(ori) * angular_vel;
  linear_vel_dot = -1 * Eigen::Vector3d(0, 0, constants::g);
  angular_vel_dot = -1 * angular_vel.cross(inertia_mat * angular_vel);
  int ridge_idx = 0;
  for(const auto & contact : motion_param.contact_list)
  {
    for(const auto & vertex_with_ridge : contact->vertexWithRidgeList_)
    {
      const auto & vertex = vertex_with_ridge.vertex;
      for(const auto & ridge : vertex_with_ridge.ridgeList)
      {
        linear_vel_dot += u[ridge_idx] * ridge / mass_;
        angular_vel_dot += u[ridge_idx] * (vertex - pos).cross(ridge);
        ridge_idx++;
      }
    }
  }
  inertia_mat.llt().solveInPlace(angular_vel_dot);

  return x + dt_ * x_dot;
}

double DdpSingleRigidBody::DdpProblem::runningCost(double t, const StateDimVector & x, const InputDimVector & u) const
{
  const auto & ref_data = ref_data_func_(t);

  return 0.5 * weight_param_.running_pos.dot((x.segment<3>(0) - ref_data.pos).cwiseAbs2())
         + 0.5 * weight_param_.running_ori.dot((x.segment<3>(3) - ref_data.ori).cwiseAbs2())
         + 0.5 * weight_param_.running_linear_vel.dot((x.segment<3>(6)).cwiseAbs2())
         + 0.5 * weight_param_.running_angular_vel.dot((x.segment<3>(9)).cwiseAbs2())
         + 0.5 * weight_param_.running_force * u.squaredNorm();
}

double DdpSingleRigidBody::DdpProblem::terminalCost(double t, const StateDimVector & x) const
{
  const auto & ref_data = ref_data_func_(t);

  return 0.5 * weight_param_.terminal_pos.dot((x.segment<3>(0) - ref_data.pos).cwiseAbs2())
         + 0.5 * weight_param_.terminal_ori.dot((x.segment<3>(3) - ref_data.ori).cwiseAbs2())
         + 0.5 * weight_param_.terminal_linear_vel.dot((x.segment<3>(6)).cwiseAbs2())
         + 0.5 * weight_param_.terminal_angular_vel.dot((x.segment<3>(9)).cwiseAbs2());
}

void DdpSingleRigidBody::DdpProblem::calcStateEqDeriv(double t,
                                                      const StateDimVector & x,
                                                      const InputDimVector & u,
                                                      Eigen::Ref<StateStateDimMatrix> state_eq_deriv_x,
                                                      Eigen::Ref<StateInputDimMatrix> state_eq_deriv_u) const
{
  const MotionParam & motion_param = motion_param_func_(t);
  const Eigen::Ref<const Eigen::Matrix3d> inertia_mat = motion_param.inertia_mat;
  Eigen::LLT<Eigen::Matrix3d> llt;
  llt.compute(inertia_mat);

  const Eigen::Ref<const Eigen::Vector3d> & pos = x.segment<3>(0);
  const Eigen::Ref<const Eigen::Vector3d> & ori = x.segment<3>(3);
  const Eigen::Ref<const Eigen::Vector3d> & angular_vel = x.segment<3>(9);

  state_eq_deriv_x.setZero();
  state_eq_deriv_u.setZero();

  state_eq_deriv_x.block<3, 3>(0, 6).diagonal().setConstant(1.0);
  state_eq_deriv_x.block<3, 3>(3, 9) = matAngularVelToEulerDot(ori);

  // The symbolic mathematics library (SymPy) was used to get the following code
  double w1 = angular_vel[0];
  double w2 = angular_vel[1];
  double w3 = angular_vel[2];
  double cos_alpha = std::cos(ori.x());
  double sin_alpha = std::sin(ori.x());
  double cos_beta = std::cos(ori.y());
  double sin_beta = std::sin(ori.y());
  double cos_beta_2 = std::pow(cos_beta, 2);
  double sin_beta_2 = std::pow(sin_beta, 2);
  double I11 = inertia_mat(0, 0);
  double I12 = inertia_mat(0, 1);
  double I13 = inertia_mat(0, 2);
  double I22 = inertia_mat(1, 1);
  double I23 = inertia_mat(1, 2);
  double I33 = inertia_mat(2, 2);
  state_eq_deriv_x.block<3, 3>(3, 3).col(0)
      << -w1 * sin_alpha * sin_beta / cos_beta + w2 * sin_beta * cos_alpha / cos_beta,
      -w1 * cos_alpha - w2 * sin_alpha, -w1 * sin_alpha / cos_beta + w2 * cos_alpha / cos_beta;
  state_eq_deriv_x.block<3, 3>(3, 3).col(1) << w1 * sin_beta_2 * cos_alpha / cos_beta_2 + w1 * cos_alpha
                                                   + w2 * sin_alpha * sin_beta_2 / cos_beta_2 + w2 * sin_alpha,
      0.0, w1 * sin_beta * cos_alpha / cos_beta_2 + w2 * sin_alpha * sin_beta / cos_beta_2;
  state_eq_deriv_x.block<3, 3>(9, 9) << I12 * w3 - I13 * w2, -I13 * w1 + I22 * w3 - 2 * I23 * w2 - I33 * w3,
      I12 * w1 + I22 * w2 + 2 * I23 * w3 - I33 * w2, -I11 * w3 + 2 * I13 * w1 + I23 * w2 + I33 * w3,
      -I12 * w3 + I23 * w1, -I11 * w1 - I12 * w2 - 2 * I13 * w3 + I33 * w1,
      I11 * w2 - 2 * I12 * w1 - I22 * w2 - I23 * w3, I11 * w1 + 2 * I12 * w2 + I13 * w3 - I22 * w1, I13 * w2 - I23 * w1;
  llt.solveInPlace(state_eq_deriv_x.block<3, 3>(9, 9));

  int ridge_idx = 0;
  Eigen::Vector3d totalForce = Eigen::Vector3d::Zero();
  for(const auto & contact : motion_param.contact_list)
  {
    for(const auto & vertex_with_ridge : contact->vertexWithRidgeList_)
    {
      const auto & vertex = vertex_with_ridge.vertex;
      for(const auto & ridge : vertex_with_ridge.ridgeList)
      {
        totalForce += u[ridge_idx] * ridge;
        state_eq_deriv_u.middleRows<3>(6).col(ridge_idx) = ridge / mass_;
        state_eq_deriv_u.middleRows<3>(9).col(ridge_idx) = (vertex - pos).cross(ridge);
        ridge_idx++;
      }
    }
  }
  state_eq_deriv_x.block<3, 3>(9, 0) = llt.solve(crossMat(totalForce));
  llt.solveInPlace(state_eq_deriv_u.middleRows<3>(9));

  state_eq_deriv_x *= dt_;
  state_eq_deriv_x.diagonal().array() += 1.0;
  state_eq_deriv_u *= dt_;
}

void DdpSingleRigidBody::DdpProblem::calcRunningCostDeriv(double t,
                                                          const StateDimVector & x,
                                                          const InputDimVector & u,
                                                          Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                                          Eigen::Ref<InputDimVector> running_cost_deriv_u) const
{
  const auto & ref_data = ref_data_func_(t);

  running_cost_deriv_x << weight_param_.running_pos.cwiseProduct(x.segment<3>(0) - ref_data.pos),
      weight_param_.running_ori.cwiseProduct(x.segment<3>(3) - ref_data.ori),
      weight_param_.running_linear_vel.cwiseProduct(x.segment<3>(6)),
      weight_param_.running_angular_vel.cwiseProduct(x.segment<3>(9));
  running_cost_deriv_u = weight_param_.running_force * u;
}

void DdpSingleRigidBody::DdpProblem::calcRunningCostDeriv(double t,
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
  running_cost_deriv_xx.diagonal() << weight_param_.running_pos, weight_param_.running_ori,
      weight_param_.running_linear_vel, weight_param_.running_angular_vel;
  running_cost_deriv_uu.setIdentity();
  running_cost_deriv_uu.diagonal() *= weight_param_.running_force;
  running_cost_deriv_xu.setZero();
}

void DdpSingleRigidBody::DdpProblem::calcTerminalCostDeriv(double t,
                                                           const StateDimVector & x,
                                                           Eigen::Ref<StateDimVector> terminal_cost_deriv_x) const
{
  const auto & ref_data = ref_data_func_(t);

  terminal_cost_deriv_x << weight_param_.terminal_pos.cwiseProduct(x.segment<3>(0) - ref_data.pos),
      weight_param_.terminal_ori.cwiseProduct(x.segment<3>(3) - ref_data.ori),
      weight_param_.terminal_linear_vel.cwiseProduct(x.segment<3>(6)),
      weight_param_.terminal_angular_vel.cwiseProduct(x.segment<3>(9));
}

void DdpSingleRigidBody::DdpProblem::calcTerminalCostDeriv(double t,
                                                           const StateDimVector & x,
                                                           Eigen::Ref<StateDimVector> terminal_cost_deriv_x,
                                                           Eigen::Ref<StateStateDimMatrix> terminal_cost_deriv_xx) const
{
  calcTerminalCostDeriv(t, x, terminal_cost_deriv_x);

  terminal_cost_deriv_xx.setZero();
  terminal_cost_deriv_xx.diagonal() << weight_param_.terminal_pos, weight_param_.terminal_ori,
      weight_param_.terminal_linear_vel, weight_param_.terminal_angular_vel;
}

DdpSingleRigidBody::InitialParam::InitialParam(const DdpProblem::StateDimVector & state)
{
  pos << state.segment<3>(0);
  ori << state.segment<3>(3);
  linear_vel << state.segment<3>(6);
  angular_vel << state.segment<3>(9);
}

DdpSingleRigidBody::DdpProblem::StateDimVector DdpSingleRigidBody::InitialParam::toState() const
{
  DdpProblem::StateDimVector state;
  state << pos, ori, linear_vel, angular_vel;
  return state;
}

DdpSingleRigidBody::DdpSingleRigidBody(double mass,
                                       double horizon_dt,
                                       int horizon_steps,
                                       const WeightParam & weight_param)
: ddp_problem_(std::make_shared<DdpProblem>(horizon_dt, mass, weight_param)),
  ddp_solver_(std::make_shared<nmpc_ddp::DDPSolver<12, Eigen::Dynamic>>(ddp_problem_))
{
  ddp_solver_->config().with_input_constraint = true;
  ddp_solver_->config().horizon_steps = horizon_steps;
  ddp_solver_->config().initial_lambda = 1e-6;
  ddp_solver_->config().lambda_min = 1e-8;
  ddp_solver_->config().lambda_thre = 1e-7;
  ddp_solver_->setInputLimitsFunc([this](double t) -> std::array<Eigen::VectorXd, 2> {
    std::array<Eigen::VectorXd, 2> limits;
    int input_dim = ddp_problem_->inputDim(t);
    limits[0].setConstant(input_dim, force_scale_limits_[0]);
    limits[1].setConstant(input_dim, force_scale_limits_[1]);
    return limits;
  });
}

Eigen::VectorXd DdpSingleRigidBody::planOnce(const std::function<MotionParam(double)> & motion_param_func,
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
    ddp_solver_->solve(current_time, initial_param.toState(), current_u_list);
  }
  else
  {
    ddp_solver_->solve(current_time, initial_param.toState(), initial_param.u_list);
  }

  return ddp_solver_->controlData().u_list[0];
}
