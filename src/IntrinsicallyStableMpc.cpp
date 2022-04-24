/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/IntrinsicallyStableMpc.h>

using namespace CCC;

IntrinsicallyStableMpc1d::IntrinsicallyStableMpc1d(double com_height,
                                                   double horizon_duration,
                                                   double horizon_dt,
                                                   QpSolverCollection::QpSolverType qp_solver_type,
                                                   const WeightParam & weight_param)
: weight_param_(weight_param), horizon_dt_(horizon_dt),
  horizon_steps_(static_cast<int>(std::ceil(horizon_duration / horizon_dt))),
  omega_(std::sqrt(constants::g / com_height)), lambda_(std::exp(-1 * omega_ * horizon_dt_))
{
  // Set P
  P_.setZero(horizon_steps_, horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    for(int j = 0; j < i + 1; j++)
    {
      P_(i, j) = horizon_dt_;
    }
  }

  // Setup QP
  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
  qp_coeff_.setup(horizon_steps_, 1, 2 * horizon_steps_);
  qp_coeff_.obj_mat_.diagonal().setConstant(weight_param_.zmp_vel);
  // Although no formula is given in the paper, it is mentioned that reference ZMP is considered in the objective
  // function.
  qp_coeff_.obj_mat_.noalias() += weight_param_.zmp * P_.transpose() * P_;
  qp_coeff_.obj_vec_.setZero();
  // See equation (14) in the paper
  qp_coeff_.eq_mat_(0, 0) = (1 - lambda_) / (omega_ * (1 - std::pow(lambda_, horizon_steps_)));
  for(int i = 1; i < horizon_steps_; i++)
  {
    qp_coeff_.eq_mat_(0, i) = lambda_ * qp_coeff_.eq_mat_(0, i - 1);
  }
  // See equation (8) in the paper
  qp_coeff_.ineq_mat_ << -1 * P_, P_;
  qp_coeff_.x_min_.setConstant(-1e10);
  qp_coeff_.x_max_.setConstant(1e10);
}

double IntrinsicallyStableMpc1d::planOnce(const std::function<RefData(double)> & ref_data_func,
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

double IntrinsicallyStableMpc1d::procOnce(const std::vector<RefData> & ref_data_seq,
                                          const InitialParam & initial_param,
                                          double current_time,
                                          double control_dt)
{
  std::array<double, 2> current_zmp_limits;

  // Set QP coefficients
  // See equation (14) in the paper
  qp_coeff_.eq_vec_ << initial_param.capture_point - initial_param.planned_zmp;
  // See equation (8) in the paper
  Eigen::VectorXd ref_zmp_vec(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    const auto & ref_data = ref_data_seq[i];
    ref_zmp_vec[i] = ref_data.zmp;
    qp_coeff_.ineq_vec_[i] = -1 * ref_data.zmp_limits[0];
    qp_coeff_.ineq_vec_[i + horizon_steps_] = ref_data.zmp_limits[1];

    if(i == 0)
    {
      current_zmp_limits = ref_data.zmp_limits;
    }
  }
  qp_coeff_.obj_vec_ = weight_param_.zmp * P_.transpose()
                       * (Eigen::VectorXd::Constant(horizon_steps_, initial_param.planned_zmp) - ref_zmp_vec);
  qp_coeff_.ineq_vec_.head(horizon_steps_).array() += initial_param.planned_zmp;
  qp_coeff_.ineq_vec_.tail(horizon_steps_).array() -= initial_param.planned_zmp;

  // Solve QP
  double zmp_vel = qp_solver_->solve(qp_coeff_)[0];

  // Calculate ZMP
  if(control_dt < 0)
  {
    control_dt = horizon_dt_;
  }
  double zmp =
      std::clamp(initial_param.planned_zmp + control_dt * zmp_vel, current_zmp_limits[0], current_zmp_limits[1]);

  return zmp;
}

Eigen::Vector2d IntrinsicallyStableMpc::planOnce(const std::function<RefData(double)> & ref_data_func,
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
      ref_data_seq_x_[i].zmp = ref_data.zmp.x();
      ref_data_seq_x_[i].zmp_limits[j] = ref_data.zmp_limits[j].x();
      ref_data_seq_y_[i].zmp = ref_data.zmp.y();
      ref_data_seq_y_[i].zmp_limits[j] = ref_data.zmp_limits[j].y();
    }
  }

  // Calculate ZMP
  Eigen::Vector2d planned_zmp;

  IntrinsicallyStableMpc1d::InitialParam initial_param_x;
  initial_param_x.capture_point = initial_param.capture_point.x();
  initial_param_x.planned_zmp = initial_param.planned_zmp.x();
  planned_zmp.x() = mpc_1d_->procOnce(ref_data_seq_x_, initial_param_x, current_time, control_dt);

  IntrinsicallyStableMpc1d::InitialParam initial_param_y;
  initial_param_y.capture_point = initial_param.capture_point.y();
  initial_param_y.planned_zmp = initial_param.planned_zmp.y();
  planned_zmp.y() = mpc_1d_->procOnce(ref_data_seq_y_, initial_param_y, current_time, control_dt);

  return planned_zmp;
}
