/* Author: Masaki Murooka */

#include <CCC/PreviewControlZmp.h>

using namespace CCC;

PreviewControl<3, 1, 1>::WeightParam PreviewControlZmp::WeightParam::toPreviewControlWeightParam() const
{
  PreviewControl<3, 1, 1>::WeightParam pc_weight_param;
  pc_weight_param.output << zmp;
  pc_weight_param.input << com_jerk;
  return pc_weight_param;
}

double PreviewControlZmp::planOnce(const std::function<double(double)> & ref_zmp_func,
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
