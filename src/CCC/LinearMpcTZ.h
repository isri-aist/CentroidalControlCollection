/* Author: Masaki Murooka */

#pragma once

#include <qp_solver_collection/QpSolver.h>

#include <CCC/EigenTypes.h>
#include <CCC/VariantSequentialExtension.h>

namespace CCC
{
/** \brief Linear MPC of translational z-component motion. */
class LinearMpcTZ
{
public:
  /** \brief State dimension. */
  static constexpr int state_dim_ = 2;

  /** \brief Type of state-space model. */
  using _StateSpaceModel = StateSpaceModel<state_dim_, Eigen::Dynamic, Eigen::Dynamic>;

  /** \brief Type of state vector. */
  using StateDimVector = Eigen::Matrix<double, state_dim_, 1>;

public:
  /** \brief State-space model for contact phase. */
  class ModelContactPhase : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    ModelContactPhase(double mass);
  };

  /** \brief State-space model for non-contact phase. */
  class ModelNoncontactPhase : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    ModelNoncontactPhase(double mass);
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param mass robot mass [kg]
      \param dt discretization timestep [sec]
      \param qp_solver_type QP solver type
  */
  LinearMpcTZ(double mass,
              double dt,
              QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::QLD);

  /** \brief Plan one step.
      \param contact_seq sequence of contact/non-contact phases (true for contact phase)
      \param current_pos current position
      \param current_vel current velocity
      \param ref_pos_seq sequence of reference position (same length as contact_seq)
      \param pos_weight objective weight for position
      \param force_weight objective weight for force
      \returns planned force sequence
      \note The length of contact_seq is assumed to be horizon size.
  */
  Eigen::VectorXd runOnce(const std::vector<bool> & contact_seq,
                          double current_pos,
                          double current_vel,
                          const Eigen::VectorXd & ref_pos_seq,
                          double pos_weight = 1.0,
                          double force_weight = 1e-7);

  /** \brief Plan with loop.
      \param contact_seq sequence of contact/non-contact phases (true for contact phase)
      \param current_pos current position
      \param current_vel current velocity
      \param ref_pos_seq sequence of reference position (same length as contact_seq)
      \param horizon_duration horizon duration (i.e., duration of preview window) [sec]
      \param pos_weight objective weight for position
      \param force_weight objective weight for force
      \returns planned force sequence
  */
  void runLoop(const std::vector<bool> & contact_seq,
               double current_pos,
               double current_vel,
               const Eigen::VectorXd & ref_pos_seq,
               double horizon_duration = 4.0,
               double pos_weight = 1.0,
               double force_weight = 1e-7);

protected:
  /** \brief Plan one step.
      \tparam ListType type of state-space model list
      \param model_list state-space model list
      \param current_x current state (product of mass and position, linear momentum)
      \param ref_pos_seq sequence of reference position
      \param pos_weight objective weight for position
      \param force_weight objective weight for force
      \returns planned force sequence
  */
  template<template<class> class ListType>
  Eigen::VectorXd runOnce(const ListType<std::shared_ptr<_StateSpaceModel>> & model_list,
                          const StateDimVector & current_x,
                          const Eigen::VectorXd & ref_pos_seq,
                          double pos_weight,
                          double force_weight);

public:
  //! Robot mass [kg]
  double mass_ = 0;

  //! Discretization timestep [sec]
  double dt_ = 0;

  //! State-space model for contact phase
  std::shared_ptr<_StateSpaceModel> model_contact_phase_;

  //! State-space model for non-contact phase
  std::shared_ptr<_StateSpaceModel> model_noncontact_phase_;

  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;

  //! Min/max z-component force [N]
  std::pair<double, double> force_range_ = {10.0, 4000.0};

  //! Planned position sequence
  Eigen::VectorXd planned_pos_seq_;

  //! Planned velocity sequence
  Eigen::VectorXd planned_vel_seq_;

  //! Planned force sequence
  Eigen::VectorXd planned_force_seq_;
};
} // namespace CCC
