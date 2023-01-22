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

  state_eq_deriv_x.setZero();
  state_eq_deriv_x.block<3, 3>(0, 3).diagonal().setConstant(1 / mass_);
  state_eq_deriv_x.block<3, 3>(6, 0) = crossMat(ridges_mat * u);
  state_eq_deriv_x *= dt_;
  state_eq_deriv_x.diagonal().array() += 1.0;

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

DdpCentroidal::InitialParam::InitialParam(const DdpProblem::StateDimVector & state, double mass)
{
  pos << state.segment<3>(0);
  vel << state.segment<3>(3) / mass;
  angular_momentum << state.segment<3>(6);
}

DdpCentroidal::DdpProblem::StateDimVector DdpCentroidal::InitialParam::toState(double mass) const
{
  DdpProblem::StateDimVector state;
  state << pos, mass * vel, angular_momentum;
  return state;
}

void DdpCentroidal::SimModel::procOnce(double t,
                                       const StateDimVector & x,
                                       const InputDimVector & u,
                                       Eigen::Ref<StateDimVector> next_x,
                                       Eigen::Ref<OutputDimVector> y) const
{
  const MotionParam & motion_param = motion_param_func_(t);
  const Eigen::Ref<const Eigen::Matrix3Xd> & vertices_mat = motion_param.vertex_ridge_list.topRows<3>();
  const Eigen::Ref<const Eigen::Matrix3Xd> & ridges_mat = motion_param.vertex_ridge_list.bottomRows<3>();

  const Eigen::Ref<const Eigen::Vector3d> & pos = x.segment<3>(0);
  const Eigen::Ref<const Eigen::Vector3d> & linear_momentum = x.segment<3>(3);
  const Eigen::Ref<const Eigen::Vector3d> & angular_momentum = x.segment<3>(6);

  // Set next_x
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
  next_x = x + dt_ * x_dot;

  // Set y
  Eigen::Ref<Eigen::Vector3d> pos_output = y.segment<3>(0);
  Eigen::Ref<Eigen::Vector3d> vel_output = y.segment<3>(3);
  Eigen::Ref<Eigen::Vector3d> acc_output = y.segment<3>(6);
  Eigen::Ref<Eigen::Vector3d> angular_momentum_output = y.segment<3>(9);
  pos_output = pos;
  vel_output = linear_momentum / mass_;
  acc_output = linear_momentum_dot / mass_;
  angular_momentum_output = angular_momentum;
}

DdpCentroidal::DdpCentroidal(double mass, double horizon_dt, int horizon_steps, const WeightParam & weight_param)
: ddp_problem_(std::make_shared<DdpProblem>(horizon_dt, mass, weight_param)),
  ddp_solver_(std::make_shared<nmpc_ddp::DDPSolver<9, Eigen::Dynamic>>(ddp_problem_))
{
  ddp_solver_->config().with_input_constraint = true;
  ddp_solver_->config().horizon_steps = horizon_steps;
  ddp_solver_->setInputLimitsFunc([this](double t) -> std::array<Eigen::VectorXd, 2> {
    std::array<Eigen::VectorXd, 2> limits;
    int input_dim = ddp_problem_->inputDim(t);
    limits[0].setConstant(input_dim, force_scale_limits_[0]);
    limits[1].setConstant(input_dim, force_scale_limits_[1]);
    return limits;
  });
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

void DdpCentroidal::planLoop(const std::function<MotionParam(double)> & motion_param_func,
                             const std::function<RefData(double)> & ref_data_func,
                             const InitialParam & initial_param,
                             const std::pair<double, double> & motion_time_range,
                             double sim_dt,
                             int ddp_max_iter)
{
  int seq_len = static_cast<int>((motion_time_range.second - motion_time_range.first) / sim_dt);

  // Loop
  double current_t = motion_time_range.first;
  DdpProblem::StateDimVector current_x = initial_param.toState(ddp_problem_->mass_);
  for(int i = 0; i < seq_len; i++)
  {
    // Set sim_model
    const auto & sim_model = std::make_shared<SimModel>(ddp_problem_->mass_, motion_param_func, sim_dt);

    // Calculate optimal force
    InitialParam current_initial_param(current_x, ddp_problem_->mass_);
    current_initial_param.u_list = ddp_solver_->controlData().u_list;
    if(!current_initial_param.u_list.empty())
    {
      for(int i = 0; i < ddp_solver_->config().horizon_steps; i++)
      {
        double t = current_t + i * ddp_problem_->dt();
        int input_dim = ddp_problem_->inputDim(t);
        if(current_initial_param.u_list[i].size() != input_dim)
        {
          current_initial_param.u_list[i].setZero(input_dim);
        }
      }
    }
    const Eigen::VectorXd & current_u = planOnce(motion_param_func, ref_data_func, current_initial_param, current_t);
    ddp_solver_->config().max_iter = ddp_max_iter; // Set max_iter from second simulation iteration

    // Save current data
    DdpProblem::StateDimVector next_x;
    SimModel::OutputDimVector current_y;
    sim_model->procOnce(current_t, current_x, current_u, next_x, current_y);
    const auto & current_ref_data = ref_data_func(current_t);
    MotionData current_motion_data;
    current_motion_data.time = current_t;
    current_motion_data.ref_pos = current_ref_data.pos;
    current_motion_data.ref_vel << 0.0, 0.0, 0.0;
    current_motion_data.ref_acc << 0.0, 0.0, 0.0;
    current_motion_data.ref_angular_momentum << 0.0, 0.0, 0.0;
    current_motion_data.planned_pos = current_y.segment<3>(0);
    current_motion_data.planned_vel = current_y.segment<3>(3);
    current_motion_data.planned_acc = current_y.segment<3>(6);
    current_motion_data.planned_force = current_u;
    current_motion_data.planned_angular_momentum = current_y.segment<3>(9);
    motion_data_seq_.emplace(current_t, current_motion_data);

    // Simulate one step
    current_t += sim_dt;
    current_x = next_x;
  }
}

void DdpCentroidal::dumpMotionDataSeq(const std::string & file_path, bool print_command) const
{
  std::ofstream ofs(file_path);
  ofs << "time ref_pos_x ref_pos_y ref_pos_z ref_vel_x ref_vel_y ref_vel_z ref_angular_momentum_x "
         "ref_angular_momentum_y ref_angular_momentum_z planned_pos_x planned_pos_y planned_pos_z planned_vel_x "
         "planned_vel_y planned_vel_z planned_acc_x planned_acc_y planned_acc_z planned_angular_momentum_x "
         "planned_angular_momentum_y planned_angular_momentum_z planned_force"
      << std::endl;
  for(const auto & motion_data_kv : motion_data_seq_)
  {
    motion_data_kv.second.dump(ofs);
  }
  if(print_command)
  {
    std::cout << "Run the following commands in gnuplot:\n"
              << "  set key autotitle columnhead\n"
              << "  set key noenhanced\n"
              << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:11 w lp # pos_x\n"
              << "  plot \"" << file_path << "\" u 1:4 w lp, \"\" u 1:13 w lp # pos_z\n"
              << "  plot \"" << file_path << "\" u 1:9 w lp, \"\" u 1:21 w lp # momentum_y\n";
  }
}
