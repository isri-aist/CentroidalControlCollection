/* Author: Masaki Murooka */

#include <algorithm>

#include <CCC/LinearMpcZmp.h>

using namespace CCC;

LinearMpcZmp1d::LinearMpcZmp1d(double com_height,
                               double horizon_duration,
                               double horizon_dt,
                               QpSolverCollection::QpSolverType qp_solver_type)
: horizon_dt_(horizon_dt), horizon_steps_(static_cast<int>(std::ceil(horizon_duration / horizon_dt))),
  model_(std::make_shared<ComZmpModelJerkInput>(com_height))
{
  // Setup model
  model_->calcDiscMatrix(horizon_dt_);
  seq_ext_ = std::make_shared<InvariantSequentialExtension<3, 1, 1>>(model_, horizon_steps_, true);

  // Setup QP
  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
  qp_coeff_.setup(horizon_steps_, 0, 2 * horizon_steps_);
  qp_coeff_.obj_mat_.setIdentity();
  qp_coeff_.obj_vec_.setZero();
  qp_coeff_.ineq_mat_ << -1 * seq_ext_->B_seq_, seq_ext_->B_seq_;
  qp_coeff_.x_min_.setConstant(-1e10);
  qp_coeff_.x_max_.setConstant(1e10);
}

double LinearMpcZmp1d::planOnce(const std::function<RefData(double)> & ref_data_func,
                                const InitialParam & initial_param,
                                double current_time,
                                double control_dt)
{
  // Set ref_data_seq
  std::vector<RefData> ref_data_seq(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    double t = current_time + i * horizon_dt_;
    ref_data_seq[i] = ref_data_func(t);
  }

  return procOnce(ref_data_seq, initial_param, current_time, control_dt);
}

double LinearMpcZmp1d::procOnce(const std::vector<RefData> & ref_data_seq,
                                const InitialParam & initial_param,
                                double, // current_time
                                double control_dt)
{
  std::array<double, 2> current_zmp_limits;

  // Set QP coefficients
  qp_coeff_.ineq_vec_.head(horizon_steps_) = seq_ext_->A_seq_ * initial_param;
  qp_coeff_.ineq_vec_.tail(horizon_steps_) = -1 * qp_coeff_.ineq_vec_.head(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    const auto & ref_data = ref_data_seq[i];
    qp_coeff_.ineq_vec_[i] -= ref_data.zmp_limits[0];
    qp_coeff_.ineq_vec_[i + horizon_steps_] += ref_data.zmp_limits[1];

    if(i == 0)
    {
      current_zmp_limits = ref_data.zmp_limits;
    }
  }

  // Solve QP
  double com_jerk = qp_solver_->solve(qp_coeff_)[0];

  // Calculate ZMP
  if(control_dt < 0)
  {
    control_dt = horizon_dt_;
  }
  double com_acc = initial_param[2] + control_dt * com_jerk;
  double com_pos = initial_param[0] + control_dt * initial_param[1] + 0.5 * std::pow(control_dt, 2) * initial_param[2];
  double zmp = std::clamp(com_pos + model_->C_(0, 2) * com_acc, current_zmp_limits[0], current_zmp_limits[1]);

  return zmp;
}

Eigen::Vector2d LinearMpcZmp::planOnce(const std::function<RefData(double)> & ref_data_func,
                                       const InitialParam & initial_param,
                                       double current_time,
                                       double control_dt)
{
  // Set ref_data_seq
  for(int i = 0; i < mpc_1d_->horizon_steps_; i++)
  {
    double t = current_time + i * mpc_1d_->horizon_dt_;
    const auto & ref_data = ref_data_func(t);
    for(int j = 0; j < 2; j++)
    {
      ref_data_seq_x_[i].zmp_limits[j] = ref_data.zmp_limits[j].x();
      ref_data_seq_y_[i].zmp_limits[j] = ref_data.zmp_limits[j].y();
    }
  }

  // Calculate ZMP
  Eigen::Vector2d planned_zmp;

  LinearMpcZmp1d::InitialParam initial_param_x;
  initial_param_x << initial_param.pos.x(), initial_param.vel.x(), initial_param.acc.x();
  planned_zmp.x() = mpc_1d_->procOnce(ref_data_seq_x_, initial_param_x, current_time, control_dt);

  LinearMpcZmp1d::InitialParam initial_param_y;
  initial_param_y << initial_param.pos.y(), initial_param.vel.y(), initial_param.acc.y();
  planned_zmp.y() = mpc_1d_->procOnce(ref_data_seq_y_, initial_param_y, current_time, control_dt);

  return planned_zmp;
}
