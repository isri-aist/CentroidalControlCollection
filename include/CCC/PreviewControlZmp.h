/* Author: Masaki Murooka */

#pragma once

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
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
  };

  /** \brief State-space model of CoM-ZMP dynamics with CoM jerk input and ZMP output. */
  class ComZmpModel : public StateSpaceModel<3, 1, 1>
  {
  public:
    /** \brief Constructor.
        \param com_height height of robot CoM [m]
    */
    ComZmpModel(double com_height)
    {
      A_(0, 1) = 1;
      A_(1, 2) = 1;

      B_(2, 0) = 1;

      C_(0, 0) = 1;
      C_(0, 2) = -1 * com_height / constants::g;
    }
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param model state-space model
   */
  PreviewControlZmp(const std::shared_ptr<ComZmpModel> & model) : PreviewControl<3, 1, 1>(model) {}

  /** \brief Calculate the gain of the preview control.
      \param horizon horizon of the preview control [sec]
      \param dt sampling time of the preview control [sec]
      \param weight_param objective weight parameter
   */
  inline void calcGain(double horizon, double dt, const WeightParam & weight_param = WeightParam())
  {
    Eigen::Vector1d output_weight;
    output_weight << weight_param.zmp;
    Eigen::Vector1d input_weight;
    input_weight << weight_param.com_jerk;
    PreviewControl::calcGain(horizon, dt, output_weight, input_weight);
  }

  /** \brief Plan ZMP.
      \param initial_param initial parameter (CoM position, velocity, and acceleration)
      \param ref_zmp_traj trajectory of reference ZMP
      \return planned ZMP
   */
  inline double planZmp(const InitialParam & initial_param, const Eigen::VectorXd & ref_zmp_traj) const
  {
    Eigen::Vector1d optimal_u = PreviewControl::calcOptimalInput(initial_param, ref_zmp_traj);
    PreviewControl<3, 1, 1>::StateDimVector optimal_x = model_->stateEqDisc(initial_param, optimal_u);
    return model_->observEq(optimal_x)[0];
  }
};
} // namespace CCC
