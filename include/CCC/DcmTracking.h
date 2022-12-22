/* Author: Masaki Murooka */

#pragma once

#include <map>

#include <Eigen/Core>

#include <CCC/Constants.h>

namespace CCC
{
/** \brief Walking control based on tracking of divergent component of motion (DCM).

    \todo Support three-dimensional motion.

    \todo Is it possible to handle continuous ZMP transition during the double support phase?

    See the following for a detailed formulation.
      - J Englsberger, et al. Three-dimensional bipedal walking control using divergent component of motion. IROS, 2013.
 */
class DcmTracking
{
public:
  /** \brief Reference data. */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Current ZMP [m]
    Eigen::Vector2d current_zmp = Eigen::Vector2d::Zero();

    /** \brief List of pairs of future ZMP and switching time

        In the referenced paper, the length of time_zmp_list is three.
    */
    std::map<double, Eigen::Vector2d> time_zmp_list;
  };

  /** \brief Initial parameter.

      Initial parameter is DCM only.
  */
  using InitialParam = Eigen::Vector2d;

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param feedback_gain feedback gain to calculate ZMP
   */
  DcmTracking(double com_height, double feedback_gain = 2.0)
  : feedback_gain_(feedback_gain), omega_(std::sqrt(constants::g / com_height))
  {
  }

  /** \brief Plan one step.
      \param ref_data reference data
      \param initial_param initial parameter (current DCM)
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned ZMP

      In the referenced paper, the length of time_zmp_list is three.
 */
  Eigen::Vector2d planOnce(const RefData & ref_data, const InitialParam & initial_param, double current_time) const;

public:
  //! Feedback gain to calculate control ZMP
  double feedback_gain_ = 0;

protected:
  //! Time constant for inverted pendulum dynamics
  double omega_ = 0;
};
} // namespace CCC
