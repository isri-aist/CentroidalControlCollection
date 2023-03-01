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

  /** \brief State-space model for simulation.

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
      y &= \begin{bmatrix} \dfrac{1}{m} & 0 \\ 0 & \dfrac{1}{m} \\ 0 & 0 \end{bmatrix} \boldsymbol{x} +
      \begin{bmatrix} 0 \\ 0 \\ \dfrac{1}{m} \end{bmatrix} \boldsymbol{u} +
      \begin{bmatrix} 0 \\ 0 \\ - g \end{bmatrix}
      \f}

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} c_z \\ \dot{c}_z \\ \ddot{c}_z \end{bmatrix},
      u = f_z, y = c_z
      \f}
   */
  class SimModel : public StateSpaceModel<state_dim_, 1, 3>
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    SimModel(double mass);
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param mass robot mass [kg]
      \param horizon_dt discretization timestep in horizon [sec]
      \param qp_solver_type QP solver type
  */
  LinearMpcZ(double mass,
             double horizon_dt,
             QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any);

  /** \brief Plan one step.
      \param contact_func function of contact/non-contact phases (returns true for contact phase)
      \param ref_pos_func function of reference position [m]
      \param initial_param initial parameter
      \param horizon_time_range start and end time of horizon ([sec], [sec])
      \param weight_param objective weight parameter
      \returns planned force sequence
  */
  Eigen::VectorXd planOnce(const std::function<bool(double)> & contact_func,
                           const std::function<double(double)> & ref_pos_func,
                           const InitialParam & initial_param,
                           const std::pair<double, double> & horizon_time_range,
                           const WeightParam & weight_param = WeightParam());

protected:
  /** \brief Process one step. */
  Eigen::VectorXd procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                           const StateDimVector & current_x,
                           const Eigen::VectorXd & ref_pos_seq,
                           const WeightParam & weight_param);

public:
  //! Robot mass [kg]
  double mass_ = 0;

  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

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
