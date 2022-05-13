/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolverCollection.h>

namespace CCC
{
/** \brief QP-based MPC with stability constraint for one-dimensional CoM-ZMP model.

    See the following for a detailed formulation.
      - N Scianca, et al. Intrinsically Stable MPC for Humanoid Gait Generation. Humanoids, 2016.
 */
class IntrinsicallyStableMpc1d
{
  friend class IntrinsicallyStableMpc;

public:
  /** \brief Reference data. */
  struct RefData
  {
    //! ZMP [m]
    double zmp = 0;

    //! Min/max limits of ZMP [m]
    std::array<double, 2> zmp_limits;
  };

  /** \brief Initial parameter. */
  struct InitialParam
  {
    //! Capture point [m]
    double capture_point = 0;

    //! Current ZMP planned in previous step [m]
    double planned_zmp = 0;
  };

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! ZMP weight
    double zmp;

    //! ZMP velocity weight
    double zmp_vel;

    /** \brief Constructor.
        \param _zmp ZMP weight
        \param _zmp_vel ZMP velocity weight
     */
    WeightParam(double _zmp = 1.0, double _zmp_vel = 1e-3) : zmp(_zmp), zmp_vel(_zmp_vel) {}
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param qp_solver_type QP solver type
      \param weight_param objective weight parameter
   */
  IntrinsicallyStableMpc1d(double com_height,
                           double horizon_duration,
                           double horizon_dt,
                           QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any,
                           const WeightParam & weight_param = WeightParam());

  /** \brief Plan one step.
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \param control_dt control timestep used to calculate ZMP (if omitted, horizon_dt is used)
      \returns planned ZMP
   */
  double planOnce(const std::function<RefData(double)> & ref_data_func,
                  const InitialParam & initial_param,
                  double current_time,
                  double control_dt = -1);

protected:
  /** \brief Process one step. */
  double procOnce(const std::vector<RefData> & ref_data_seq,
                  const InitialParam & initial_param,
                  double current_time,
                  double control_dt);

protected:
  //! Weight parameter
  WeightParam weight_param_;

  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

  //! Number of steps in horizon
  int horizon_steps_ = -1;

  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;

  //! Time constant for inverted pendulum dynamics
  double omega_ = 0;

  //! Constant value dependent on omega and horizon_dt
  double lambda_ = 0;

  //! Constant matrix dependent on horizon_dt (defined in equation (7) in the paper)
  Eigen::MatrixXd P_;
};

/** \brief QP-based MPC with stability constraint for CoM-ZMP model.

    See the following for a detailed formulation.
      - N Scianca, et al. Intrinsically Stable MPC for Humanoid Gait Generation. Humanoids, 2016.
 */
class IntrinsicallyStableMpc
{
public:
  /** \brief Reference data.

      \todo It is assumed that ZMP limits are independent for the x and y components. This assumption is not valid if
     the foot is placed diagonally during the single-support phase or if the feet are not aligned during the
     double-support phase.

      \todo Footstep online planning is not supported.
   */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! ZMP [m]
    Eigen::Vector2d zmp = Eigen::Vector2d::Zero();

    //! Min/max limits of ZMP [m]
    std::array<Eigen::Vector2d, 2> zmp_limits;
  };

  /** \brief Initial parameter. */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Capture point [m]
    Eigen::Vector2d capture_point = Eigen::Vector2d::Zero();

    //! Current ZMP planned in previous step [m]
    Eigen::Vector2d planned_zmp = Eigen::Vector2d::Zero();
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param qp_solver_type QP solver type
      \param weight_param objective weight parameter
   */
  IntrinsicallyStableMpc(
      double com_height,
      double horizon_duration,
      double horizon_dt,
      QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any,
      const IntrinsicallyStableMpc1d::WeightParam & weight_param = IntrinsicallyStableMpc1d::WeightParam())
  : mpc_1d_(
      std::make_shared<IntrinsicallyStableMpc1d>(com_height, horizon_duration, horizon_dt, qp_solver_type, weight_param))
  {
    ref_data_seq_x_.resize(mpc_1d_->horizon_steps_);
    ref_data_seq_y_.resize(mpc_1d_->horizon_steps_);
  }

  /** \brief Plan one step.
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \param control_dt control timestep used to calculate ZMP (if omitted, horizon_dt is used)
      \returns planned ZMP
   */
  Eigen::Vector2d planOnce(const std::function<RefData(double)> & ref_data_func,
                           const InitialParam & initial_param,
                           double current_time,
                           double control_dt = -1);

protected:
  //! One-dimensional linear MPC
  std::shared_ptr<IntrinsicallyStableMpc1d> mpc_1d_;

  //! Reference data sequence of x
  std::vector<IntrinsicallyStableMpc1d::RefData> ref_data_seq_x_;

  //! Reference data sequence of y
  std::vector<IntrinsicallyStableMpc1d::RefData> ref_data_seq_y_;
};
} // namespace CCC
