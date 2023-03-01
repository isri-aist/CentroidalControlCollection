/* Author: Masaki Murooka */

#include <fstream>

#include <ForceColl/Contact.h>

#include <CCC/Constants.h>
#include <CCC/LinearMpcXY.h>

using namespace CCC;

namespace
{
/** \brief Get the number of ridges from contact list. */
int getTotalRidgeNum(const std::vector<std::shared_ptr<ForceColl::Contact>> & contact_list)
{
  int ridge_num = 0;
  for(const auto & contact : contact_list)
  {
    ridge_num += contact->ridgeNum();
  }
  return ridge_num;
}
} // namespace

LinearMpcXY::StateDimVector LinearMpcXY::InitialParam::toState(double mass) const
{
  StateDimVector state;
  state << mass * pos.x(), mass * vel.x(), mass * pos.y(), mass * vel.y(), angular_momentum;
  return state;
}

LinearMpcXY::StateDimVector LinearMpcXY::RefData::toOutput(double mass) const
{
  StateDimVector output;
  output << mass * pos.x(), mass * vel.x(), mass * pos.y(), mass * vel.y(), angular_momentum;
  return output;
}

Eigen::VectorXd LinearMpcXY::WeightParam::inputWeight(int total_input_dim) const
{
  return Eigen::VectorXd::Constant(total_input_dim, force);
}

Eigen::VectorXd LinearMpcXY::WeightParam::outputWeight(size_t seq_len) const
{
  StateDimVector output_weight_one;
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
: StateSpaceModel(LinearMpcXY::state_dim_, getTotalRidgeNum(motion_param.contact_list), output_dim),
  motion_param_(motion_param)
{
  A_(0, 1) = 1;
  A_(2, 3) = 1;
  A_(4, 2) = -1 * motion_param_.total_force_z / mass;
  A_(5, 0) = motion_param_.total_force_z / mass;

  int ridge_idx = 0;
  for(const auto & contact : motion_param_.contact_list)
  {
    for(const auto & vertex_with_ridge : contact->vertexWithRidgeList_)
    {
      const auto & vertex = vertex_with_ridge.vertex;
      for(const auto & ridge : vertex_with_ridge.ridgeList)
      {
        B_.col(ridge_idx) << 0, ridge.x(), 0, ridge.y(),
            -1 * (vertex.z() - motion_param_.com_z) * ridge.y() + vertex.y() * ridge.z(),
            (vertex.z() - motion_param_.com_z) * ridge.x() + -1 * vertex.x() * ridge.z();
        ridge_idx++;
      }
    }
  }
}

LinearMpcXY::LinearMpcXY(double mass,
                         double horizon_dt,
                         int horizon_steps,
                         const WeightParam & weight_param,
                         QpSolverCollection::QpSolverType qp_solver_type)
: mass_(mass), horizon_dt_(horizon_dt), horizon_steps_(horizon_steps), weight_param_(weight_param),
  force_range_(3.0, 3.0 * mass * constants::g)
{
  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

Eigen::VectorXd LinearMpcXY::planOnce(const std::function<MotionParam(double)> & motion_param_func,
                                      const std::function<RefData(double)> & ref_data_func,
                                      const InitialParam & initial_param,
                                      double current_time)
{
  // Set model_list and ref_output_seq
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list(horizon_steps_);
  Eigen::VectorXd ref_output_seq(horizon_steps_ * RefData::outputDim());
  for(int i = 0; i < horizon_steps_; i++)
  {
    double t = current_time + i * horizon_dt_;
    model_list[i] = std::make_shared<Model>(mass_, motion_param_func(t));
    model_list[i]->calcDiscMatrix(horizon_dt_);
    ref_output_seq.segment<RefData::outputDim()>(i * RefData::outputDim()) = ref_data_func(t).toOutput(mass_);
  }

  // Calculate optimal force
  return procOnce(model_list, initial_param.toState(mass_), ref_output_seq);
}

Eigen::VectorXd LinearMpcXY::procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                                      const StateDimVector & current_x,
                                      const Eigen::VectorXd & ref_output_seq)
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
  const Eigen::VectorXd & output_weight = weight_param_.outputWeight(model_list.size());
  qp_coeff_.obj_mat_.noalias() = seq_ext.B_seq_.transpose() * output_weight.asDiagonal() * seq_ext.B_seq_;
  // Diagonal matrix cannot be added at the same time.
  // See https://forum.kde.org/viewtopic.php?f=74&t=136617#p365547
  qp_coeff_.obj_mat_.diagonal() += weight_param_.inputWeight(total_input_dim);
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

    int ridge_idx = 0;
    for(const auto & contact : model->motion_param_.contact_list)
    {
      for(const auto & vertex_with_ridge : contact->vertexWithRidgeList_)
      {
        for(const auto & ridge : vertex_with_ridge.ridgeList)
        {
          qp_coeff_.eq_mat_(accum_eq_dim, accum_input_dim + ridge_idx) = ridge.z();
          ridge_idx++;
        }
      }
    }
    qp_coeff_.eq_vec_[accum_eq_dim] = model->motion_param_.total_force_z;

    accum_eq_dim++;
    accum_input_dim += model->inputDim();
  }
  qp_coeff_.x_min_.setConstant(force_range_.first);
  qp_coeff_.x_max_.setConstant(force_range_.second);

  // Solve QP
  return qp_solver_->solve(qp_coeff_).head(model_list[0]->inputDim());
}
