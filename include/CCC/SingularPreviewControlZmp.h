/* Author: Masaki Murooka */

#pragma once

#include <functional>
#include <memory>

#include <Eigen/Dense>

#include <CCC/Constants.h>

namespace CCC
{
/** \brief Singular preview control for one-dimensional CoM-ZMP model.

    See the following for a detailed formulation.
      - J Urata, et al. Online Decision of Foot Placement using Singular LQ Preview Regulation. Humanoids, 2011.
 */
class SingularPreviewControlZmp1d
{
  friend class SingularPreviewControlZmp;

public:
  /** \brief Initial parameter. */
  struct InitialParam
  {
    //! CoM position [m]
    double pos = 0;

    //! CoM velocity [m/s]
    double vel = 0;

    //! Current ZMP planned in previous step [m]
    double planned_zmp = 0;
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
   */
  SingularPreviewControlZmp1d(double com_height, double horizon_duration, double horizon_dt)
  : horizon_dt_(horizon_dt), horizon_steps_(static_cast<int>(std::ceil(horizon_duration / horizon_dt))),
    omega_(std::sqrt(constants::g / com_height))
  {
  }

  /** \brief Plan one step.
      \param ref_zmp_func function of reference ZMP [m]
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \param control_dt control timestep used to calculate ZMP (if omitted, horizon_dt is used)
      \returns planned ZMP
   */
  double planOnce(const std::function<double(double)> & ref_zmp_func,
                  const InitialParam & initial_param,
                  double current_time,
                  double control_dt = -1) const;

protected:
  /** \brief Process one step. */
  double procOnce(const Eigen::VectorXd & ref_zmp_seq,
                  const InitialParam & initial_param,
                  double current_time,
                  double control_dt) const;

protected:
  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

  //! Number of steps in horizon
  int horizon_steps_ = -1;

  //! Time constant for inverted pendulum dynamics
  double omega_ = 0;
};

/** \brief Singular preview control for CoM-ZMP model.

    See the following for a detailed formulation.
      - J Urata, et al. Online Decision of Foot Placement using Singular LQ Preview Regulation. Humanoids, 2011.
 */
class SingularPreviewControlZmp
{
public:
  /** \brief Initial parameter. */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector2d pos = Eigen::Vector2d::Zero();

    //! CoM velocity [m/s]
    Eigen::Vector2d vel = Eigen::Vector2d::Zero();

    //! Current ZMP planned in previous step [m]
    Eigen::Vector2d planned_zmp = Eigen::Vector2d::Zero();
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
   */
  SingularPreviewControlZmp(double com_height, double horizon_duration, double horizon_dt)
  : spc_1d_(std::make_shared<SingularPreviewControlZmp1d>(com_height, horizon_duration, horizon_dt))
  {
  }

  /** \brief Plan one step.
      \param ref_zmp_func function of reference ZMP [m]
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \param control_dt control timestep used to calculate ZMP (if omitted, horizon_dt is used)
      \returns planned ZMP
   */
  Eigen::Vector2d planOnce(const std::function<Eigen::Vector2d(double)> & ref_zmp_func,
                           const InitialParam & initial_param,
                           double current_time,
                           double control_dt = -1) const;

public:
  //! One-dimensional preview control
  std::shared_ptr<SingularPreviewControlZmp1d> spc_1d_;
};
} // namespace CCC
