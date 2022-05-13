/* Author: Masaki Murooka */

#include <CCC/SingularPreviewControlZmp.h>

using namespace CCC;

double SingularPreviewControlZmp1d::planOnce(const std::function<double(double)> & ref_zmp_func,
                                             const InitialParam & initial_param,
                                             double current_time,
                                             double control_dt) const
{
  // Set ref_zmp_seq
  Eigen::VectorXd ref_zmp_seq(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    double t = current_time + i * horizon_dt_;
    ref_zmp_seq[i] = ref_zmp_func(t);
  }

  return procOnce(ref_zmp_seq, initial_param, current_time, control_dt);
}

double SingularPreviewControlZmp1d::procOnce(const Eigen::VectorXd & ref_zmp_seq,
                                             const InitialParam & initial_param,
                                             double current_time,
                                             double control_dt) const
{
  // Calculate feedback term
  Eigen::Vector3d x;
  x << initial_param.planned_zmp, initial_param.pos, initial_param.vel;
  // Equation (16) in the paper
  Eigen::Vector3d K;
  K << (1 + omega_ * horizon_dt_) / horizon_dt_ + omega_ / (1 + omega_ * horizon_dt_),
      -1 * (2 + omega_ * horizon_dt_) / horizon_dt_, -1 * (2 + omega_ * horizon_dt_) / (omega_ * horizon_dt_);
  double u_fb = -1 * K.dot(x);

  // Calculate feedforward term
  // Equation (27) in the paper
  double S0 = ref_zmp_seq[ref_zmp_seq.size() - 1] * (1 + omega_ * horizon_dt_) / (omega_ * horizon_dt_);
  for(int i = ref_zmp_seq.size() - 2; i >= 1; i--)
  {
    // Equation (25) in the paper
    S0 = ref_zmp_seq[i] + S0 / (1 + omega_ * horizon_dt_);
  }
  // Equation (26) in the paper
  double u_ff =
      ref_zmp_seq[0] / horizon_dt_ - omega_ * (2 + omega_ * horizon_dt_) * S0 / std::pow(1 + omega_ * horizon_dt_, 2);

  // Calculate ZMP velocity
  // Equation (9) in the paper
  double u = u_fb + u_ff;

  // Calculate ZMP
  double zmp = initial_param.planned_zmp + control_dt * u;

  return zmp;
}

Eigen::Vector2d SingularPreviewControlZmp::planOnce(const std::function<Eigen::Vector2d(double)> & ref_zmp_func,
                                                    const InitialParam & initial_param,
                                                    double current_time,
                                                    double control_dt) const
{
  // Set ref_zmp_seq
  Eigen::Matrix2Xd ref_zmp_seq(2, spc_1d_->horizon_steps_);
  for(int i = 0; i < spc_1d_->horizon_steps_; i++)
  {
    double t = current_time + i * spc_1d_->horizon_dt_;
    ref_zmp_seq.col(i) = ref_zmp_func(t);
  }

  // Calculate ZMP
  Eigen::Vector2d planned_zmp;

  SingularPreviewControlZmp1d::InitialParam initial_param_x;
  initial_param_x.pos = initial_param.pos.x();
  initial_param_x.vel = initial_param.vel.x();
  initial_param_x.planned_zmp = initial_param.planned_zmp.x();
  planned_zmp.x() = spc_1d_->procOnce(ref_zmp_seq.row(0), initial_param_x, current_time, control_dt);

  SingularPreviewControlZmp1d::InitialParam initial_param_y;
  initial_param_y.pos = initial_param.pos.y();
  initial_param_y.vel = initial_param.vel.y();
  initial_param_y.planned_zmp = initial_param.planned_zmp.y();
  planned_zmp.y() = spc_1d_->procOnce(ref_zmp_seq.row(1), initial_param_y, current_time, control_dt);

  return planned_zmp;
}
