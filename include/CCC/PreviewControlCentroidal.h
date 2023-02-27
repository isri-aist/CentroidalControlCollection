/* Author: Masaki Murooka */

#pragma once

#include <mc_rtc/Configuration.h>

#include <CCC/EigenTypes.h>
#include <CCC/PreviewControl.h>

namespace ForceColl
{
class Contact;
}

namespace CCC
{
/** \brief Preview control for one-dimensional centroidal model.

    See the following for a detailed formulation.
      - M Murooka, et al. Centroidal trajectory generation and stabilization based on preview control for humanoid
   multi-contact motion. RA-Letters, 2022.
 */
class PreviewControlCentroidal1d : public PreviewControl<3, 1, 2>
{
  friend class PreviewControlCentroidal;

public:
  /** \brief Reference data.

      For the position components, the first element is the CoM position and the second element is the force; for the
     orientation components, the first element is the base link orientation and the second element is the moment.
  */
  using RefData = Eigen::Vector2d;

  /** \brief Initial parameter.

      For the position components, it consists of the position, velocity, and acceleration of the CoM. For the
     orientation components, it consists of the orientation, angular velocity, and angular acceleration of the base
     link.
  */
  using InitialParam = PreviewControl<3, 1, 2>::StateDimVector;

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! Position weight
    double pos;

    //! Wrench weight
    double wrench;

    //! Jerk weight
    double jerk;

    /** \brief Constructor.
        \param _pos position weight
        \param _wrench wrench weight
        \param _jerk jerk weight
     */
    WeightParam(double _pos = 2e2, double _wrench = 5e-4, double _jerk = 1e-8) : pos(_pos), wrench(_wrench), jerk(_jerk)
    {
    }

    /** \brief Convert to PreviewControl::WeightParam. */
    PreviewControl<3, 1, 2>::WeightParam toPreviewControlWeightParam() const;
  };

  /** \brief State-space model of one-dimensional centroidal dynamics with jerk input and position and wrench output.

      Dynamics is expressed by the following equation.
      \f{align*}{
      m \ddot{c} = f
      \f}
      \f$c\f$ is position (CoM position or base link orientation).
      \f$f\f$ is wrench (force for the position components, moment for the orientation components).
      \f$m\f$ is inertia parameter (mass for the position components, moment of inertia for the orientation components).

      This can be represented as a linear time-invariant system as follows.
      \f{align*}{
      \boldsymbol{\dot{x}} &=
      \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix}
      \boldsymbol{x} +
      \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} u \\
      \boldsymbol{y} &= \begin{bmatrix} 1 & 0 & 0 \\ 0 & 0 & m \end{bmatrix} \boldsymbol{x}
      \f}

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} c \\ \dot{c} \\ \ddot{c} \end{bmatrix},
      u = \dddot{c}, \boldsymbol{y} = \begin{bmatrix} c \\ f \end{bmatrix}
      \f}
   */
  class CentroidalModel1d : public StateSpaceModel<3, 1, 2>
  {
  public:
    /** \brief Constructor.
        \param inertia_param inertia parameter (mass [kg] for the position components, moment of inertia [kg m^2] for
       the orientation components)
    */
    CentroidalModel1d(double inertia_param);
  };

public:
  /** \brief Constructor.
      \param inertia_param inertia parameter (mass [kg] for the position components, moment of inertia [kg m^2] for the
             orientation components)
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControlCentroidal1d(double inertia_param,
                             double horizon_duration,
                             double horizon_dt,
                             const WeightParam & weight_param = WeightParam())
  : PreviewControl<3, 1, 2>(std::make_shared<CentroidalModel1d>(inertia_param),
                            horizon_duration,
                            horizon_dt,
                            weight_param.toPreviewControlWeightParam())
  {
  }

protected:
  /** \brief Process one step. */
  double procOnce(const Eigen::VectorXd & ref_output_seq,
                  const InitialParam & initial_param,
                  double current_time,
                  double control_dt) const;
};

/** \brief Preview control for centroidal model.

    See the following for a detailed formulation.
      - M Murooka, et al. Centroidal trajectory generation and stabilization based on preview control for humanoid
   multi-contact motion. RA-Letters, 2022.
 */
class PreviewControlCentroidal
{
public:
  /** \brief Motion parameter. */
  struct MotionParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** \brief Contact list. */
    std::vector<std::shared_ptr<ForceColl::Contact>> contact_list;
  };

  /** \brief Reference data.

      The first three elements are position components and the latter three are orientation components.
   */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m] and base link orientation [rad]
    Eigen::Vector6d pos = Eigen::Vector6d::Zero();

    //! Force [N] and moment [Nm]
    Eigen::Vector6d wrench = Eigen::Vector6d::Zero();
  };

  /** \brief Initial parameter.

      The first three elements are position components and the latter three are orientation components.
   */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m] and base link orientation [rad]
    Eigen::Vector6d pos = Eigen::Vector6d::Zero();

    //! CoM velocity [m/s] and base link angular velocity [rad/s]
    Eigen::Vector6d vel = Eigen::Vector6d::Zero();

    //! CoM acceleration [m/s^2] and base link angular acceleration [rad/s^2]
    Eigen::Vector6d acc = Eigen::Vector6d::Zero();

    /** \brief Convert to PreviewControlCentroidal1d::InitialParam.
        \param idx component index
     */
    PreviewControlCentroidal1d::InitialParam toInitialParam1d(int idx) const;
  };

  /** \brief Weight parameter. */
  struct WeightParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Position weight
    Eigen::Vector6d pos;

    //! Wrench weight
    Eigen::Vector6d wrench;

    //! Jerk weight
    Eigen::Vector6d jerk;

    //! Configuration of wrench distribution
    mc_rtc::Configuration wrench_dist_config;

    /** \brief Constructor.
        \param _pos position weight
        \param _wrench wrench weight
        \param _jerk jerk weight
        \param _wrench_dist_config configuration of wrench distribution
     */
    WeightParam(const Eigen::Vector6d & _pos = (Eigen::Vector6d() << 2e2, 2e2, 2e2, 1e2, 1e2, 1e2).finished(),
                const Eigen::Vector6d & _wrench = (Eigen::Vector6d() << 5e-4, 5e-4, 5e-4, 5e-3, 5e-3, 5e-3).finished(),
                const Eigen::Vector6d & _jerk = Eigen::Vector6d::Constant(1e-8),
                const mc_rtc::Configuration & _wrench_dist_config = {})
    : pos(_pos), wrench(_wrench), jerk(_jerk), wrench_dist_config(_wrench_dist_config)
    {
    }

    /** \brief Convert to PreviewControlCentroidal1d::WeightParam.
        \param idx component index
     */
    PreviewControlCentroidal1d::WeightParam toWeightParam1d(int idx) const;
  };

public:
  /** \brief Constructor.
      \param mass mass [kg]
      \param moment_of_inertia moment of inertia [kg m^2]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControlCentroidal(double mass,
                           const Eigen::Vector3d & moment_of_inertia,
                           double horizon_duration,
                           double horizon_dt,
                           const WeightParam & weight_param = WeightParam());

  /** \brief Plan one step.
      \param motion_param motion parameter
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \param control_dt control timestep used to calculate the planned wrench (if omitted, horizon_dt is used)
      \returns planned wrench
   */
  Eigen::Vector6d planOnce(const MotionParam & motion_param,
                           const std::function<RefData(double)> & ref_data_func,
                           const InitialParam & initial_param,
                           double current_time,
                           double control_dt = -1) const;

public:
  //! One-dimensional preview control
  std::array<std::shared_ptr<PreviewControlCentroidal1d>, 6> preview_control_1d_;

  //! Configuration of wrench distribution
  mc_rtc::Configuration wrench_dist_config_;

protected:
  //! Mass [kg]
  double mass_;
};
} // namespace CCC
