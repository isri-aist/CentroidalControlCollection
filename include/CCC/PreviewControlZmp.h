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

  /** \brief State-space model of CoM-ZMP dynamics with CoM jerk input and ZMP output.

      Dynamics is expressed by the following equation.
      \f{align*}{
      \ddot{c}_x = \dfrac{g}{c_z} (c_x - z_x)
      \f}
      \f$\boldsymbol{c}\f$ and \f$\boldsymbol{z}\f$ are CoM and ZMP.

      This can be represented as a linear time-invariant system as follows.
      \f{align*}{
      \boldsymbol{\dot{x}} &=
      \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix}
      \boldsymbol{x} +
      \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} u \\
      y &= \begin{bmatrix} 1 & 0 & - \dfrac{c_z}{g} \end{bmatrix} \boldsymbol{x}
      \f}

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} c_x \\ \dot{c}_x \\ \ddot{c}_x \end{bmatrix},
      u = \dddot{c}_x, y = z_x
      \f}
   */
  class Model : public StateSpaceModel<3, 1, 1>
  {
  public:
    /** \brief Constructor.
        \param com_height height of robot CoM [m]
    */
    Model(double com_height);
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControlZmp(double com_height,
                    double horizon_duration,
                    double horizon_dt,
                    const WeightParam & weight_param = WeightParam())
  : PreviewControl<3, 1, 1>(std::make_shared<Model>(com_height),
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
};
} // namespace CCC
