/* Author: Masaki Murooka */

#include <fstream>

#include <CCC/Constants.h>
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

LinearMpcZ::LinearMpcZ(double mass,
                       double horizon_dt,
                       int horizon_steps,
                       const WeightParam & weight_param,
                       QpSolverCollection::QpSolverType qp_solver_type)
: mass_(mass), horizon_dt_(horizon_dt), horizon_steps_(horizon_steps), weight_param_(weight_param),
  force_range_(10.0, 10.0 * mass * constants::g)
{
  model_contact_ = std::make_shared<ModelContactPhase>(mass_);
  model_noncontact_ = std::make_shared<ModelNoncontactPhase>(mass_);

  model_contact_->calcDiscMatrix(horizon_dt_);
  model_noncontact_->calcDiscMatrix(horizon_dt_);

  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

double LinearMpcZ::planOnce(const std::function<bool(double)> & contact_func,
                            const std::function<double(double)> & ref_pos_func,
                            const InitialParam & initial_param,
                            double current_time)
{
  // Planned force is always zero if there is no contact
  if(!contact_func(current_time))
  {
    return 0.0;
  }

  // Set model_list and ref_pos_seq
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list(horizon_steps_);
  Eigen::VectorXd ref_pos_seq(horizon_steps_);
  for(int i = 0; i < horizon_steps_; i++)
  {
    double t = current_time + i * horizon_dt_;
    model_list[i] = contact_func(t) ? model_contact_ : model_noncontact_;
    ref_pos_seq[i] = ref_pos_func(t);
  }

  // Calculate optimal force
  return procOnce(model_list, mass_ * initial_param, ref_pos_seq);
}

double LinearMpcZ::procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                            const StateDimVector & current_x,
                            const Eigen::VectorXd & ref_pos_seq)
{
  // Calculate sequential extension
  VariantSequentialExtension<state_dim_> seq_ext(model_list, true);

  // Set QP coefficients
  int total_input_dim = seq_ext.totalInputDim();
  if(!(qp_coeff_.dim_var_ == total_input_dim && qp_coeff_.dim_eq_ == 0 && qp_coeff_.dim_ineq_ == 0))
  {
    qp_coeff_.setup(total_input_dim, 0, 0);
  }
  qp_coeff_.obj_mat_.noalias() = weight_param_.pos * seq_ext.B_seq_.transpose() * seq_ext.B_seq_;
  qp_coeff_.obj_mat_.diagonal().array() += weight_param_.force;
  qp_coeff_.obj_vec_.noalias() =
      -1 * weight_param_.pos * seq_ext.B_seq_.transpose() * (ref_pos_seq - seq_ext.A_seq_ * current_x - seq_ext.E_seq_);
  qp_coeff_.x_min_.setConstant(force_range_.first);
  qp_coeff_.x_max_.setConstant(force_range_.second);

  // Solve QP
  return qp_solver_->solve(qp_coeff_)[0];
}
