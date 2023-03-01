/* Author: Masaki Murooka */

#include <ForceColl/WrenchDistribution.h>

#include <CCC/Constants.h>
#include <CCC/PreviewControlCentroidal.h>

using namespace CCC;

PreviewControlCentroidal1d::CentroidalModel1d::CentroidalModel1d(double inertia_param)
{
  A_(0, 1) = 1;
  A_(1, 2) = 1;

  B_(2, 0) = 1;

  C_(0, 0) = 1;
  C_(1, 2) = inertia_param;
}

PreviewControl<3, 1, 2>::WeightParam PreviewControlCentroidal1d::WeightParam::toPreviewControlWeightParam() const
{
  PreviewControl<3, 1, 2>::WeightParam pc_weight_param;
  pc_weight_param.output << pos, wrench;
  pc_weight_param.input << jerk;
  return pc_weight_param;
}

double PreviewControlCentroidal1d::procOnce(const Eigen::VectorXd & ref_output_seq,
                                            const InitialParam & initial_param,
                                            double, // current_time
                                            double control_dt) const
{
  // Calculate jerk
  double jerk = PreviewControl::calcOptimalInput(initial_param, ref_output_seq)[0];

  // Calculate next state and output
  if(control_dt < 0)
  {
    control_dt = horizon_dt_;
  }
  double acc = initial_param[2] + control_dt * jerk;
  double force = model_->C_(1, 2) * acc;

  return force;
}

PreviewControlCentroidal1d::InitialParam PreviewControlCentroidal::InitialParam::toInitialParam1d(int idx) const
{
  return PreviewControlCentroidal1d::InitialParam(pos[idx], vel[idx], acc[idx]);
}

PreviewControlCentroidal1d::WeightParam PreviewControlCentroidal::WeightParam::toWeightParam1d(int idx) const
{
  return PreviewControlCentroidal1d::WeightParam(pos[idx], wrench[idx], jerk[idx]);
}

PreviewControlCentroidal::PreviewControlCentroidal(double mass,
                                                   const Eigen::Vector3d & moment_of_inertia,
                                                   double horizon_duration,
                                                   double horizon_dt,
                                                   const PreviewControlCentroidal::WeightParam & weight_param)
: mass_(mass)
{
  for(int i = 0; i < 6; i++)
  {
    double inertia_param = (i < 3 ? mass_ : moment_of_inertia[i - 3]);
    preview_control_1d_[i] = std::make_shared<PreviewControlCentroidal1d>(inertia_param, horizon_duration, horizon_dt,
                                                                          weight_param.toWeightParam1d(i));
  }

  wrench_dist_config_ = weight_param.wrench_dist_config;
}

Eigen::Vector6d PreviewControlCentroidal::planOnce(const MotionParam & motion_param,
                                                   const std::function<RefData(double)> & ref_data_func,
                                                   const InitialParam & initial_param,
                                                   double current_time,
                                                   double control_dt) const
{
  int horizon_steps = preview_control_1d_[0]->horizon_steps_;
  double horizon_dt = preview_control_1d_[0]->horizon_dt_;

  // Set ref_output_seq
  Eigen::Matrix<double, 6, Eigen::Dynamic> ref_output_seq(6, horizon_steps * 2);
  for(int i = 0; i < horizon_steps; i++)
  {
    double t = current_time + (i + 1) * horizon_dt;
    const auto & ref_data = ref_data_func(t);
    ref_output_seq.col(2 * i) = ref_data.pos;
    ref_output_seq.col(2 * i + 1) = ref_data.wrench;
  }

  // Run preview control
  Eigen::Vector6d planned_wrench;
  for(int i = 0; i < 6; i++)
  {
    planned_wrench[i] = preview_control_1d_[i]->procOnce(ref_output_seq.row(i).transpose(),
                                                         initial_param.toInitialParam1d(i), current_time, control_dt);
  }

  // Project wrench
  sva::ForceVecd wrench_without_gravity = sva::ForceVecd(planned_wrench.tail<3>(), planned_wrench.head<3>());
  wrench_without_gravity.force() += Eigen::Vector3d(0, 0, mass_ * constants::g);
  auto wrench_dist = std::make_shared<ForceColl::WrenchDistribution>(motion_param.contact_list, wrench_dist_config_);
  wrench_dist->run(wrench_without_gravity, initial_param.pos.head<3>());

  return (Eigen::Vector6d() << wrench_dist->resultTotalWrench_.force(), wrench_dist->resultTotalWrench_.moment())
      .finished();
}
