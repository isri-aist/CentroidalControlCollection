/* Author: Masaki Murooka */

#include <CCC/PreviewControlZmp.h>

using namespace CCC;

PreviewControl<3, 1, 1>::WeightParam PreviewControlZmp1d::WeightParam::toPreviewControlWeightParam() const
{
  PreviewControl<3, 1, 1>::WeightParam pc_weight_param;
  pc_weight_param.output << zmp;
  pc_weight_param.input << com_jerk;
  return pc_weight_param;
}

double PreviewControlZmp1d::planOnce(const std::function<double(double)> & ref_zmp_func,
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

double PreviewControlZmp1d::procOnce(const Eigen::VectorXd & ref_zmp_seq,
                                     const InitialParam & initial_param,
                                     double current_time,
                                     double control_dt) const
{
  // Calculate CoM jerk
  double com_jerk = PreviewControl::calcOptimalInput(initial_param, ref_zmp_seq)[0];

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

Eigen::Vector2d PreviewControlZmp::planOnce(const std::function<Eigen::Vector2d(double)> & ref_zmp_func,
                                            const InitialParam & initial_param,
                                            double current_time,
                                            double control_dt) const
{
  // Set ref_zmp_seq
  Eigen::Matrix2Xd ref_zmp_seq(2, preview_control_1d_->horizon_steps_);
  for(int i = 0; i < preview_control_1d_->horizon_steps_; i++)
  {
    double t = current_time + i * preview_control_1d_->horizon_dt_;
    ref_zmp_seq.col(i) = ref_zmp_func(t);
  }

  // Calculate ZMP
  Eigen::Vector2d planned_zmp;

  PreviewControlZmp1d::InitialParam initial_param_x;
  initial_param_x << initial_param.pos.x(), initial_param.vel.x(), initial_param.acc.x();
  planned_zmp.x() = preview_control_1d_->procOnce(ref_zmp_seq.row(0), initial_param_x, current_time, control_dt);

  PreviewControlZmp1d::InitialParam initial_param_y;
  initial_param_y << initial_param.pos.y(), initial_param.vel.y(), initial_param.acc.y();
  planned_zmp.y() = preview_control_1d_->procOnce(ref_zmp_seq.row(1), initial_param_y, current_time, control_dt);

  return planned_zmp;
}
