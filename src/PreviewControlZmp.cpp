/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/PreviewControlZmp.h>

using namespace CCC;

PreviewControl<3, 1, 1>::WeightParam PreviewControlZmp::WeightParam::toPreviewControlWeightParam() const
{
  PreviewControl<3, 1, 1>::WeightParam pc_weight_param;
  pc_weight_param.output << zmp;
  pc_weight_param.input << com_jerk;
  return pc_weight_param;
}

PreviewControlZmp::Model::Model(double com_height)
{
  A_(0, 1) = 1;
  A_(1, 2) = 1;

  B_(2, 0) = 1;

  C_(0, 0) = 1;
  C_(0, 2) = -1 * com_height / constants::g;
}

double PreviewControlZmp::planOnce(const std::function<double(double)> & ref_zmp_func,
                                   const InitialParam & initial_param,
                                   double horizon_start_time) const
{
  // Set ref_zmp_seq
  Eigen::VectorXd ref_zmp_seq(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    double t = horizon_start_time + i * horizon_dt_;
    ref_zmp_seq[i] = ref_zmp_func(t);
  }

  // Calculate planned ZMP
  Eigen::Vector1d optimal_u = PreviewControl::calcOptimalInput(initial_param, ref_zmp_seq);
  PreviewControl<3, 1, 1>::StateDimVector optimal_x = model_->stateEqDisc(initial_param, optimal_u);
  return model_->observEq(optimal_x)[0];
}
