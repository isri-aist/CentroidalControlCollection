/* Author: Masaki Murooka */

#include <fstream>

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/LinearMpcZ.h>

using namespace CCC;

LinearMpcZ::ModelContactPhase::ModelContactPhase(double mass) : StateSpaceModel(LinearMpcZ::state_dim_, 1, 1)
{
  A_ << 0, 1, 0, 0;

  B_ << 0, 1;

  C_ << 1 / mass, 0;

  E_ << 0, -1 * mass * constants::g;
}

LinearMpcZ::ModelNoncontactPhase::ModelNoncontactPhase(double mass) : StateSpaceModel(LinearMpcZ::state_dim_, 0, 1)
{
  A_ << 0, 1, 0, 0;

  C_ << 1 / mass, 0;

  E_ << 0, -1 * mass * constants::g;
}

LinearMpcZ::SimModel::SimModel(double mass)
{
  A_ << 0, 1, 0, 0;

  B_ << 0, 1;

  C_ << 1 / mass, 0, 0, 1 / mass, 0, 0;

  D_ << 0, 0, 1 / mass;

  E_ << 0, -1 * mass * constants::g;

  F_ << 0, 0, -1 * constants::g;
}

LinearMpcZ::LinearMpcZ(double mass, double horizon_dt, QpSolverCollection::QpSolverType qp_solver_type)
: mass_(mass), horizon_dt_(horizon_dt), force_range_(10.0, 10.0 * mass * constants::g)
{
  model_contact_ = std::make_shared<ModelContactPhase>(mass_);
  model_noncontact_ = std::make_shared<ModelNoncontactPhase>(mass_);

  model_contact_->calcDiscMatrix(horizon_dt_);
  model_noncontact_->calcDiscMatrix(horizon_dt_);

  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

Eigen::VectorXd LinearMpcZ::planOnce(const std::function<bool(double)> & contact_func,
                                     const std::function<double(double)> & ref_pos_func,
                                     const InitialParam & initial_param,
                                     const std::pair<double, double> & horizon_time_range,
                                     const WeightParam & weight_param)
{
  // Set model_list and ref_pos_seq
  int horizon_steps = static_cast<int>((horizon_time_range.second - horizon_time_range.first) / horizon_dt_);
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list(horizon_steps);
  Eigen::VectorXd ref_pos_seq(horizon_steps);
  for(int i = 0; i < horizon_steps; i++)
  {
    double t = horizon_time_range.first + i * horizon_dt_;
    model_list[i] = contact_func(t) ? model_contact_ : model_noncontact_;
    ref_pos_seq[i] = ref_pos_func(t);
  }

  // Calculate optimal force
  return procOnce(model_list, mass_ * initial_param, ref_pos_seq, weight_param);
}

Eigen::VectorXd LinearMpcZ::procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                                     const StateDimVector & current_x,
                                     const Eigen::VectorXd & ref_pos_seq,
                                     const WeightParam & weight_param)
{
  // Calculate sequential extension
  VariantSequentialExtension<state_dim_> seq_ext(model_list, true);

  // Set QP coefficients
  int total_input_dim = seq_ext.totalInputDim();
  if(!(qp_coeff_.dim_var_ == total_input_dim && qp_coeff_.dim_eq_ == 0 && qp_coeff_.dim_ineq_ == 0))
  {
    qp_coeff_.setup(total_input_dim, 0, 0);
  }
  qp_coeff_.obj_mat_.noalias() = weight_param.pos * seq_ext.B_seq_.transpose() * seq_ext.B_seq_;
  qp_coeff_.obj_mat_.diagonal().array() += weight_param.force;
  qp_coeff_.obj_vec_.noalias() =
      -1 * weight_param.pos * seq_ext.B_seq_.transpose() * (ref_pos_seq - seq_ext.A_seq_ * current_x - seq_ext.E_seq_);
  qp_coeff_.x_min_.setConstant(force_range_.first);
  qp_coeff_.x_max_.setConstant(force_range_.second);

  // Solve QP
  return qp_solver_->solve(qp_coeff_);
}
