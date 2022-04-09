/* Author: Masaki Murooka */

#pragma once

#include <CCC/Constants.h>

namespace CCC
{
/** \brief Foot-guided control.

    The method is based on the analytical solution of an optimization problem consisting of a ZMP objective and a
   capture point constraint.

    See the following for a detailed formulation.
      - T Sugihara, et al. Foot-guided agile control of a biped robot through ZMP manipulation. IROS, 2017.
      - Y Kojio, et al. Unified balance control for biped robots including modification of footsteps with angular
   momentum and falling detection based on capturability. IROS, 2019.
 */
class FootGuidedControl
{
public:
  /** \brief Reference data. */
  struct RefData
  {
    //! Transition start ZMP [m]
    double transit_start_zmp = 0;

    //! Transition end ZMP [m]
    double transit_end_zmp = 0;

    //! Transition start time [s]
    double transit_start_time = 0;

    //! Transition duration [s]
    double transit_duration = 0;
  };

  /** \brief Initial parameter.

      Initial parameter is capture point only.
  */
  using InitialParam = double;

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
   */
  FootGuidedControl(double com_height) : omega_(std::sqrt(CCC::constants::g / com_height)) {}

  /** \brief Plan one step.
      \param ref_data reference data
      \param initial_param initial parameter (capture point)
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned ZMP
   */
  double planOnce(const RefData & ref_data, const InitialParam & initial_param, double current_time) const;

protected:
  //! Time constant for inverted pendulum dynamics
  double omega_ = 0;
};
} // namespace CCC
