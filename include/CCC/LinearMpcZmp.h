/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolverCollection.h>

#include <CCC/CommonModels.h>
#include <CCC/InvariantSequentialExtension.h>

namespace CCC
{
/** \brief QP-based MPC for CoM-ZMP model.

    See the following for a detailed formulation.
      - PB Wieber. Trajectory Free Linear Model Predictive Control for Stable Walking in the Presence of Strong
   Perturbations. Humanoids, 2006.
 */
class LinearMpcZmp
{
public:
  /** \brief Reference data. */
  struct RefData
  {
    //! Min/max limits of ZMP [m]
    std::array<double, 2> zmp_limits;
  };

  /** \brief Initial parameter.

      First element is CoM position, second element is CoM velocity, and third element is CoM acceleration.
  */
  using InitialParam = Eigen::Vector3d;

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param qp_solver_type QP solver type
   */
  LinearMpcZmp(double com_height,
               double horizon_duration,
               double horizon_dt,
               QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any);

  /** \brief Plan one step.
      \param ref_data_func function of reference data
      \param initial_param initial parameter (CoM position, velocity, and acceleration)
      \param current_time current time (i.e., start time of horizon) [sec]
      \param control_dt control timestep used to calculate ZMP (if omitted, horizon_dt is used)
      \returns planned ZMP
   */
  double planOnce(const std::function<RefData(double)> & ref_data_func,
                  const InitialParam & initial_param,
                  double current_time,
                  double control_dt = -1);

protected:
  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

  //! Number of steps in horizon
  int horizon_steps_ = -1;

  //! State-space model
  std::shared_ptr<ComZmpModelJerkInput> model_;

  //! Sequential extension of state-space model
  std::shared_ptr<InvariantSequentialExtension<3, 1, 1>> seq_ext_;

  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;
};
} // namespace CCC
