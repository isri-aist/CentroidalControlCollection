/* Author: Masaki Murooka */

#include <CCC/LinearMpcZmp.h>

using namespace CCC;

LinearMpcZmp::LinearMpcZmp(double com_height,
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

double LinearMpcZmp::planOnce(const std::function<RefData(double)> & ref_data_func,
                              const InitialParam & initial_param,
                              double current_time,
                              double control_dt)
{
  // Set QP coefficients
  qp_coeff_.ineq_vec_.head(horizon_steps_) = seq_ext_->A_seq_ * initial_param;
  qp_coeff_.ineq_vec_.tail(horizon_steps_) = -1 * qp_coeff_.ineq_vec_.head(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    double t = current_time + i * horizon_dt_;
    const auto & ref_data = ref_data_func(t);
    qp_coeff_.ineq_vec_[i] -= ref_data.zmp_limits[0];
    qp_coeff_.ineq_vec_[i + horizon_steps_] += ref_data.zmp_limits[1];
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
  double zmp = com_pos + model_->C_(0, 2) * com_acc;

  return zmp;
}
