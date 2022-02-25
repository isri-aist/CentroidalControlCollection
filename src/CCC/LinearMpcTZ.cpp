/* Author: Masaki Murooka */

#include <deque>

#include <CCC/Constants.h>
#include <CCC/LinearMpcTZ.h>

using namespace CCC;

LinearMpcTZ::ModelContactPhase::ModelContactPhase(double mass) : StateSpaceModel(LinearMpcTZ::state_dim_, 1, 1)
{
  A_ << 0, 1, 0, 0;

  B_ << 0, 1;

  C_ << 1 / mass, 0;

  E_ << 0, -1 * mass * constants::g;
}

LinearMpcTZ::ModelNoncontactPhase::ModelNoncontactPhase(double mass) : StateSpaceModel(LinearMpcTZ::state_dim_, 0, 1)
{
  A_ << 0, 1, 0, 0;

  C_ << 1 / mass, 0;

  E_ << 0, -1 * mass * constants::g;
}

LinearMpcTZ::LinearMpcTZ(double mass, double dt, QpSolverCollection::QpSolverType qp_solver_type) : mass_(mass), dt_(dt)
{
  model_contact_phase_ = std::make_shared<ModelContactPhase>(mass_);
  model_noncontact_phase_ = std::make_shared<ModelNoncontactPhase>(mass_);

  model_contact_phase_->calcDiscMatrix(dt_);
  model_noncontact_phase_->calcDiscMatrix(dt_);

  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

Eigen::VectorXd LinearMpcTZ::runOnce(const std::vector<bool> & contact_seq,
                                     double current_pos,
                                     double current_vel,
                                     const Eigen::VectorXd & ref_pos_seq,
                                     double pos_weight,
                                     double force_weight)
{
  // Set model list
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list(contact_seq.size());
  for(size_t i = 0; i < contact_seq.size(); i++)
  {
    model_list[i] = contact_seq[i] ? model_contact_phase_ : model_noncontact_phase_;
  }

  StateDimVector current_x = StateDimVector(mass_ * current_pos, mass_ * current_vel);
  return runOnce(model_list, current_x, ref_pos_seq, pos_weight, force_weight);
}

void LinearMpcTZ::runLoop(const std::vector<bool> & contact_seq,
                          double current_pos,
                          double current_vel,
                          const Eigen::VectorXd & ref_pos_seq,
                          double horizon_duration,
                          double pos_weight,
                          double force_weight)
{
  // Setup model list
  int seq_len = contact_seq.size();
  int horizon_size = static_cast<int>(horizon_duration / dt_);
  std::deque<std::shared_ptr<_StateSpaceModel>> model_list;
  for(int i = 0; i < horizon_size; i++)
  {
    model_list.push_back(contact_seq[i] ? model_contact_phase_ : model_noncontact_phase_);
  }

  // Loop
  StateDimVector current_x = StateDimVector(mass_ * current_pos, mass_ * current_vel);
  Eigen::VectorXd ref_pos_seq_horizon(horizon_size);
  Eigen::VectorXd opt_force_seq(horizon_size);
  planned_pos_seq_.resize(seq_len);
  planned_vel_seq_.resize(seq_len);
  planned_force_seq_.setZero(seq_len);
  for(int i = 0; i < seq_len; i++)
  {
    const auto & current_model = model_list[0];

    // Set ref_pos_seq_horizon
    ref_pos_seq_horizon.head(std::min(horizon_size, seq_len - i)) =
        ref_pos_seq.segment(i, std::min(horizon_size, seq_len - i));
    ref_pos_seq_horizon.tail(std::max(0, i + horizon_size - seq_len)).setConstant(ref_pos_seq[seq_len - 1]);

    // Save current state
    planned_pos_seq_[i] = current_x[0] / mass_;
    planned_vel_seq_[i] = current_x[1] / mass_;

    // Calculate and save optimal force
    opt_force_seq = runOnce(model_list, current_x, ref_pos_seq_horizon, pos_weight, force_weight);
    if(current_model->inputDim() > 0)
    {
      planned_force_seq_[i] = opt_force_seq[0];
    }

    // Calculate next state
    current_x = current_model->stateEqDisc(current_x, opt_force_seq.head(current_model->inputDim()));

    // Update model_list
    model_list.pop_front();
    model_list.push_back(contact_seq[std::min(i + horizon_size, seq_len - 1)] ? model_contact_phase_
                                                                              : model_noncontact_phase_);
  }
}

template<template<class> class ListType>
Eigen::VectorXd LinearMpcTZ::runOnce(const ListType<std::shared_ptr<_StateSpaceModel>> & model_list,
                                     const StateDimVector & current_x,
                                     const Eigen::VectorXd & ref_pos_seq,
                                     double pos_weight,
                                     double force_weight)
{
  // Calculate sequential extension
  VariantSequentialExtension<state_dim_, ListType> seq_ext(model_list, true);

  // Set QP coefficients
  int total_input_dim = seq_ext.totalInputDim();
  if(!(qp_coeff_.dim_var_ == total_input_dim && qp_coeff_.dim_eq_ == 0 && qp_coeff_.dim_ineq_ == 0))
  {
    qp_coeff_.setup(total_input_dim, 0, 0);
  }
  qp_coeff_.obj_mat_.noalias() = pos_weight * seq_ext.B_seq_.transpose() * seq_ext.B_seq_
                                 + force_weight * Eigen::MatrixXd::Identity(total_input_dim, total_input_dim);
  qp_coeff_.obj_vec_.noalias() =
      -1 * pos_weight * seq_ext.B_seq_.transpose() * (ref_pos_seq - seq_ext.A_seq_ * current_x - seq_ext.E_seq_);
  qp_coeff_.x_min_.setConstant(force_range_.first);
  qp_coeff_.x_max_.setConstant(force_range_.second);

  // Solve QP
  return qp_solver_->solve(qp_coeff_);
}
