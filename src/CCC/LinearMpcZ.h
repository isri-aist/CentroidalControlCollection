/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolver.h>

#include <CCC/VariantSequentialExtension.h>

namespace CCC
{
/** \brief Linear MPC of translational z-component motion. */
class LinearMpcZ
{
public:
  /** \brief State dimension. */
  static constexpr int state_dim_ = 2;

  /** \brief Type of state-space model. */
  using _StateSpaceModel = StateSpaceModel<state_dim_, Eigen::Dynamic, Eigen::Dynamic>;

  /** \brief Type of state vector. */
  using StateDimVector = Eigen::Matrix<double, state_dim_, 1>;

public:
  /** \brief Initial parameter.

      First element is CoM position, and second element is CoM velocity.
  */
  using InitialParam = Eigen::Vector2d;

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! Position weight
    double pos;

    //! Force weight
    double force;

    /** \brief Constructor.
        \param _pos position weight
        \param _force force weight
     */
    WeightParam(double _pos = 1.0, double _force = 1e-7) : pos(_pos), force(_force) {}
  };

  /** \brief Motion data. */
  struct MotionData
  {
    //! Time [s]
    double time = 0;

    //! Contact/non-contact phase (true for contact phase)
    bool contact = 0;

    //! Reference position
    double ref_pos = 0;

    //! Planned position
    double planned_pos = 0;

    //! Planned velocity
    double planned_vel = 0;

    //! Planned acceleration
    double planned_acc = 0;

    //! Planned force
    double planned_force = 0;

    /** \brief Dump data.
        \tparam StreamType stream type
     */
    template<class StreamType>
    void dump(StreamType & ofs) const
    {
      ofs << time << " " << contact << " " << ref_pos << " " << planned_pos << " " << planned_vel << " " << planned_acc
          << " " << planned_force << std::endl;
    }
  };

  /** \brief State-space model for contact phase. */
  class ModelContactPhase : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    ModelContactPhase(double mass);
  };

  /** \brief State-space model for non-contact phase. */
  class ModelNoncontactPhase : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    ModelNoncontactPhase(double mass);
  };

  /** \brief State-space model for simulation. */
  class SimModel : public StateSpaceModel<state_dim_, 1, 3>
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
    */
    SimModel(double mass);
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param mass robot mass [kg]
      \param horizon_dt discretization timestep in horizon [s]
      \param qp_solver_type QP solver type
  */
  LinearMpcZ(double mass,
             double horizon_dt,
             QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::QLD);

  /** \brief Plan one step.
      \param contact_func function of contact/non-contact phases (returns true for contact phase)
      \param ref_pos_func function of reference position [m]
      \param initial_param initial parameter
      \param horizon_time_range start and end time of horizon ([s], [s])
      \param weight_param objective weight parameter
      \returns planned force sequence
  */
  Eigen::VectorXd planOnce(const std::function<bool(double)> & contact_func,
                           const std::function<double(double)> & ref_pos_func,
                           const InitialParam & initial_param,
                           const std::pair<double, double> & horizon_time_range,
                           const WeightParam & weight_param = WeightParam());

  /** \brief Plan with loop.
      \param contact_func function of contact/non-contact phases (returns true for contact phase)
      \param ref_pos_func function of reference position [m]
      \param initial_param initial parameter
      \param motion_time_range start and end time of motion ([s], [s])
      \param horizon_duration horizon duration [s]
      \param sim_dt discretization timestep for simulation [s]
      \param weight_param objective weight parameter
      \returns planned force sequence
  */
  void planLoop(const std::function<bool(double)> & contact_func,
                const std::function<double(double)> & ref_pos_func,
                const InitialParam & initial_param,
                const std::pair<double, double> & motion_time_range,
                double horizon_duration,
                double sim_dt,
                const WeightParam & weight_param = WeightParam());

  /** \brief Dump motion data sequence by planLoop().
      \param file_path output file path
      \param print_command whether to print the plot commands
   */
  void dumpMotionDataSeq(const std::string & file_path, bool print_command = true) const;

protected:
  /** \brief Process one step. */
  Eigen::VectorXd procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                           const StateDimVector & current_x,
                           const Eigen::VectorXd & ref_pos_seq,
                           const WeightParam & weight_param);

public:
  //! Robot mass [kg]
  double mass_ = 0;

  //! Discretization timestep in horizon [s]
  double horizon_dt_ = 0;

  //! State-space model for contact phase
  std::shared_ptr<_StateSpaceModel> model_contact_;

  //! State-space model for non-contact phase
  std::shared_ptr<_StateSpaceModel> model_noncontact_;

  //! State-space model for simulation
  std::shared_ptr<SimModel> sim_model_;

  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;

  //! Min/max z-component force [N]
  std::pair<double, double> force_range_;

  //! Motion data sequence
  std::vector<MotionData> motion_data_seq_;
};
} // namespace CCC
