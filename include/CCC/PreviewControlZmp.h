/* Author: Masaki Murooka */

#pragma once

#include <CCC/PreviewControl.h>

namespace CCC
{
/** \brief Preview control for CoM-ZMP model.

    See the following for a detailed formulation.
      - Shuuji Kajita, et al. Biped walking pattern generation by using preview control of zero-moment point. ICRA,
   2003.
 */
class PreviewControlZmp : public PreviewControl<3, 1, 1>
{
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

  /** \brief State-space model of CoM-ZMP dynamics with CoM jerk input and ZMP output. */
  class ComZmpModel : public StateSpaceModel<3, 1, 1>
  {
  public:
    /** \brief Constructor.
        \param com_height height of robot CoM [m]
    */
    ComZmpModel(double com_height);
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param model state-space model
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControlZmp(const std::shared_ptr<ComZmpModel> & model,
                    double horizon_duration,
                    double horizon_dt,
                    const WeightParam & weight_param = WeightParam())
  : PreviewControl<3, 1, 1>(model, horizon_duration, horizon_dt, weight_param.toPreviewControlWeightParam())
  {
  }

  /** \brief Plan one step.
      \param ref_zmp_func function of reference ZMP [m]
      \param initial_param initial parameter (CoM position, velocity, and acceleration)
      \param horizon_start_time start time of horizon [sec]
      \returns planned ZMP
   */
  double planOnce(const std::function<double(double)> & ref_zmp_func,
                  const InitialParam & initial_param,
                  const double & horizon_start_time) const;
};
} // namespace CCC
