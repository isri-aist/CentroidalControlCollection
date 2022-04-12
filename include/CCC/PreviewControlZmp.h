/* Author: Masaki Murooka */

#pragma once

#include <CCC/CommonModels.h>
#include <CCC/PreviewControl.h>

namespace CCC
{
/** \brief Preview control for one-dimensional CoM-ZMP model.

    See the following for a detailed formulation.
      - Shuuji Kajita, et al. Biped walking pattern generation by using preview control of zero-moment point. ICRA,
   2003.
 */
class PreviewControlZmp1d : public PreviewControl<3, 1, 1>
{
  friend class PreviewControlZmp;

public:
  /** \brief Initial parameter.

      First element is CoM position, second element is CoM velocity, and third element is CoM acceleration.
  */
  using InitialParam = PreviewControl<3, 1, 1>::StateDimVector;

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! ZMP weight
    double zmp;

    //! CoM jerk weight
    double com_jerk;

    /** \brief Constructor.
        \param _zmp ZMP weight
        \param _com_jerk CoM jerk weight
     */
    WeightParam(double _zmp = 1.0, double _com_jerk = 1e-8) : zmp(_zmp), com_jerk(_com_jerk) {}

    /** \brief Convert to PreviewControl::WeightParam. */
    PreviewControl<3, 1, 1>::WeightParam toPreviewControlWeightParam() const;
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControlZmp1d(double com_height,
                      double horizon_duration,
                      double horizon_dt,
                      const WeightParam & weight_param = WeightParam())
  : PreviewControl<3, 1, 1>(std::make_shared<ComZmpModelJerkInput>(com_height),
                            horizon_duration,
                            horizon_dt,
                            weight_param.toPreviewControlWeightParam())
  {
  }

  /** \brief Plan one step.
      \param ref_zmp_func function of reference ZMP [m]
      \param initial_param initial parameter (CoM position, velocity, and acceleration)
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
};

/** \brief Preview control for CoM-ZMP model.

    See the following for a detailed formulation.
      - Shuuji Kajita, et al. Biped walking pattern generation by using preview control of zero-moment point. ICRA,
   2003.
 */
class PreviewControlZmp
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

    //! CoM acceleration [m/s^2]
    Eigen::Vector2d acc = Eigen::Vector2d::Zero();
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControlZmp(double com_height,
                    double horizon_duration,
                    double horizon_dt,
                    const PreviewControlZmp1d::WeightParam & weight_param = PreviewControlZmp1d::WeightParam())
  : preview_control_1d_(std::make_shared<PreviewControlZmp1d>(com_height, horizon_duration, horizon_dt, weight_param))
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
  std::shared_ptr<PreviewControlZmp1d> preview_control_1d_;
};
} // namespace CCC
