/* Author: Masaki Murooka */

#include <fstream>

#include <CCC/Constants.h>
#include <CCC/LinearMpcXY.h>

using namespace CCC;

Eigen::Vector6d LinearMpcXY::InitialParam::toState(double mass) const
{
  Eigen::Vector6d state;
  state << mass * pos.x(), mass * vel.x(), mass * pos.y(), mass * vel.y(), angular_momentum;
  return state;
}

Eigen::Vector6d LinearMpcXY::RefData::toOutput(double mass) const
{
  Eigen::Vector6d output;
  output << mass * pos.x(), mass * vel.x(), mass * pos.y(), mass * vel.y(), angular_momentum;
  return output;
}

Eigen::VectorXd LinearMpcXY::WeightParam::inputWeight(int total_input_dim) const
{
  return Eigen::VectorXd::Constant(total_input_dim, force);
}

Eigen::VectorXd LinearMpcXY::WeightParam::outputWeight(size_t seq_len) const
{
  Eigen::Vector6d output_weight_one;
  output_weight_one << linear_momentum_integral.x(), linear_momentum.x(), linear_momentum_integral.y(),
      linear_momentum.y(), angular_momentum.x(), angular_momentum.y();

  Eigen::VectorXd output_weight(RefData::outputDim() * seq_len);
  for(size_t i = 0; i < seq_len; i++)
  {
    output_weight.segment<RefData::outputDim()>(i * RefData::outputDim()) << output_weight_one;
  }
  return output_weight;
}

LinearMpcXY::Model::Model(double mass, const MotionParam & motion_param, int output_dim)
: StateSpaceModel(LinearMpcXY::state_dim_, static_cast<int>(motion_param.vertex_ridge_list.cols()), output_dim),
  motion_param_(motion_param)
{
  A_(0, 1) = 1;
  A_(2, 3) = 1;
  A_(4, 2) = -1 * motion_param_.total_force_z / mass;
  A_(5, 0) = motion_param_.total_force_z / mass;

  for(int i = 0; i < motion_param_.vertex_ridge_list.cols(); i++)
  {
    const Eigen::Ref<const Eigen::Vector3d> & vertex = motion_param_.vertex_ridge_list.col(i).head<3>();
    const Eigen::Ref<const Eigen::Vector3d> & ridge = motion_param_.vertex_ridge_list.col(i).tail<3>();
    B_.col(i) << 0, ridge.x(), 0, ridge.y(),
        -1 * (vertex.z() - motion_param_.com_z) * ridge.y() + vertex.y() * ridge.z(),
        (vertex.z() - motion_param_.com_z) * ridge.x() + -1 * vertex.x() * ridge.z();
  }
}

LinearMpcXY::LinearMpcXY(double mass, double horizon_dt, QpSolverCollection::QpSolverType qp_solver_type)
: mass_(mass), horizon_dt_(horizon_dt), force_range_(3.0, 3.0 * mass * constants::g)
{
  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

Eigen::VectorXd LinearMpcXY::planOnce(const std::function<MotionParam(double)> & motion_param_func,
                                      const std::function<RefData(double)> & ref_data_func,
                                      const InitialParam & initial_param,
                                      const std::pair<double, double> & horizon_time_range,
                                      const WeightParam & weight_param)
{
  // Set model_list and ref_output_seq
  int horizon_steps = static_cast<int>((horizon_time_range.second - horizon_time_range.first) / horizon_dt_);
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list(horizon_steps);
  Eigen::VectorXd ref_output_seq(horizon_steps * RefData::outputDim());
  for(int i = 0; i < horizon_steps; i++)
  {
    double t = horizon_time_range.first + i * horizon_dt_;
    model_list[i] = std::make_shared<Model>(mass_, motion_param_func(t));
    model_list[i]->calcDiscMatrix(horizon_dt_);
    ref_output_seq.segment<RefData::outputDim()>(i * RefData::outputDim()) = ref_data_func(t).toOutput(mass_);
  }

  // Calculate optimal force
  return procOnce(model_list, initial_param.toState(mass_), ref_output_seq, weight_param);
}

Eigen::VectorXd LinearMpcXY::procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                                      const StateDimVector & current_x,
                                      const Eigen::VectorXd & ref_output_seq,
                                      const WeightParam & weight_param)
{
  // Calculate sequential extension
  VariantSequentialExtension<state_dim_> seq_ext(model_list, false);

  // Setup QP coefficients
  int total_input_dim = seq_ext.totalInputDim();
  int dim_eq = 0;
  for(const auto & model : model_list)
  {
    // Do not impose the total_force_z constraint on models with no contact
    if(model->inputDim() > 0)
    {
      dim_eq++;
    }
  }
  if(!(qp_coeff_.dim_var_ == total_input_dim && qp_coeff_.dim_eq_ == dim_eq && qp_coeff_.dim_ineq_ == 0))
  {
    qp_coeff_.setup(total_input_dim, dim_eq, 0);
  }

  // Set QP objective coefficients
  const Eigen::VectorXd & output_weight = weight_param.outputWeight(model_list.size());
  qp_coeff_.obj_mat_.noalias() = seq_ext.B_seq_.transpose() * output_weight.asDiagonal() * seq_ext.B_seq_;
  // Diagonal matrix cannot be added at the same time.
  // See https://forum.kde.org/viewtopic.php?f=74&t=136617#p365547
  qp_coeff_.obj_mat_.diagonal() += weight_param.inputWeight(total_input_dim);
  qp_coeff_.obj_vec_.noalias() = -1 * seq_ext.B_seq_.transpose() * output_weight.asDiagonal()
                                 * (ref_output_seq - seq_ext.A_seq_ * current_x - seq_ext.E_seq_);

  // Set QP constraint coefficients
  qp_coeff_.eq_mat_.setZero();
  int accum_eq_dim = 0;
  int accum_input_dim = 0;
  for(const auto & _model : model_list)
  {
    const auto & model = std::dynamic_pointer_cast<Model>(_model);
    if(model->inputDim() == 0)
    {
      continue;
    }

    for(int i = 0; i < model->inputDim(); i++)
    {
      const Eigen::Ref<const Eigen::Vector3d> & ridge = model->motion_param_.vertex_ridge_list.col(i).tail<3>();
      qp_coeff_.eq_mat_(accum_eq_dim, accum_input_dim + i) = ridge.z();
    }
    qp_coeff_.eq_vec_[accum_eq_dim] = model->motion_param_.total_force_z;

    accum_eq_dim++;
    accum_input_dim += model->inputDim();
  }
  qp_coeff_.x_min_.setConstant(force_range_.first);
  qp_coeff_.x_max_.setConstant(force_range_.second);

  // Solve QP
  return qp_solver_->solve(qp_coeff_);
}
