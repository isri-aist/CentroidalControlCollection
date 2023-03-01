/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolverCollection.h>

#include <CCC/VariantSequentialExtension.h>

namespace CCC
{
/** \brief Linear MPC of translational z-component motion. */
class LinearMpcZ
{
public:
  /** \brief State dimension. */
  static constexpr int state_dim_ = 2;

  /** \brief Type of state-space model. */
  using _StateSpaceModel = StateSpaceModel<state_dim_, Eigen::Dynamic, Eigen::Dynamic>;

  /** \brief Type of state vector. */
  using StateDimVector = _StateSpaceModel::StateDimVector;

public:
  /** \brief Initial parameter.

      First element is CoM position, and second element is CoM velocity.
  */
  using InitialParam = Eigen::Vector2d;

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! Position weight
    double pos;

    //! Force weight
    double force;

    /** \brief Constructor.
        \param _pos position weight
        \param _force force weight
     */
    WeightParam(double _pos = 1.0, double _force = 1e-7) : pos(_pos), force(_force) {}
  };

  /** \brief State-space model for contact phase.

      Dynamics is expressed by the following equation.
      \f{align*}{
      P_z &= m \dot{c}_z \\
      \dot{P}_z &= f_z - m g
      \f}
      \f$c_z\f$, \f$P_z\f$, and \f$f_z\f$ are CoM height, vertical linear momentum, and vertical contact force,
     respectively.

      This can be represented as a linear time-invariant system as follows.
      \f{align*}{
      \boldsymbol{\dot{x}} &=
      \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix}
      \boldsymbol{x} +
      \begin{bmatrix} 0 \\ 1 \end{bmatrix} u +
      \begin{bmatrix} 0 \\ - m g \end{bmatrix} \\
      y &= \begin{bmatrix} \dfrac{1}{m} & 0 \end{bmatrix} \boldsymbol{x}
      \f}

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} m c_z \\ P_z \end{bmatrix},
      u = f_z, y = c_z
      \f}
   */
  class ModelContactPhase : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    ModelContactPhase(double mass);
  };

  /** \brief State-space model for non-contact phase.

      Dynamics is expressed by the following equation.
      \f{align*}{
      P_z &= m \dot{c}_z \\
      \dot{P}_z &= - m g
      \f}
      \f$c_z\f$ and \f$P_z\f$ are CoM height and vertical linear momentum, respectively.

      This can be represented as a linear time-invariant system as follows.
      \f{align*}{
      \boldsymbol{\dot{x}} &=
      \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix}
      \boldsymbol{x} +
      \begin{bmatrix} 0 \\ - m g \end{bmatrix} \\
      y &= \begin{bmatrix} \dfrac{1}{m} & 0 \end{bmatrix} \boldsymbol{x}
      \f}

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} m c_z \\ P_z \end{bmatrix},
      u = \emptyset, y = c_z
      \f}
   */
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
      \param horizon_dt discretization timestep in horizon [sec]
      \param horizon_steps number of steps in horizon
      \param weight_param objective weight parameter
      \param qp_solver_type QP solver type
  */
  LinearMpcZ(double mass,
             double horizon_dt,
             int horizon_steps,
             const WeightParam & weight_param = WeightParam(),
             QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any);

  /** \brief Plan one step.
      \param contact_func function of contact/non-contact phases (returns true for contact phase)
      \param ref_pos_func function of reference position [m]
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned force
  */
  double planOnce(const std::function<bool(double)> & contact_func,
                  const std::function<double(double)> & ref_pos_func,
                  const InitialParam & initial_param,
                  double current_time);

protected:
  /** \brief Process one step. */
  double procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                  const StateDimVector & current_x,
                  const Eigen::VectorXd & ref_pos_seq);

public:
  //! Robot mass [kg]
  double mass_ = 0;

  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

  //! Number of steps in horizon
  int horizon_steps_ = 0;

  //! Objective weight parameter
  WeightParam weight_param_;

  //! State-space model for contact phase
  std::shared_ptr<_StateSpaceModel> model_contact_;

  //! State-space model for non-contact phase
  std::shared_ptr<_StateSpaceModel> model_noncontact_;

  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;

  //! Min/max z-component force [N]
  std::pair<double, double> force_range_;
};
} // namespace CCC
