/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolverCollection.h>

#include <CCC/EigenTypes.h>
#include <CCC/MotionData.h>
#include <CCC/VariantSequentialExtension.h>

namespace CCC
{
/** \brief Linear MPC of translational and rotational x,y-component motion. */
class LinearMpcXY
{
public:
  /** \brief State dimension. */
  static constexpr int state_dim_ = 6;

  /** \brief Type of state-space model. */
  using _StateSpaceModel = StateSpaceModel<state_dim_, Eigen::Dynamic, Eigen::Dynamic>;

  /** \brief Type of state vector. */
  using StateDimVector = Eigen::Matrix<double, state_dim_, 1>;

public:
  /** \brief Motion parameter. */
  struct MotionParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM z position [m]
    double com_z;

    //! Total z force [N]
    double total_force_z;

    //! List of contact vertex and force direction (i.e., friction pyramid ridge)
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> vertex_ridge_list;
  };

  /** \brief Initial parameter. */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector2d pos = Eigen::Vector2d::Zero();

    //! CoM linear velocity [m/s]
    Eigen::Vector2d vel = Eigen::Vector2d::Zero();

    //! Angular momentum [kg m^2/s]
    Eigen::Vector2d angular_momentum = Eigen::Vector2d::Zero();

    /** \brief Get state of state-space model.
        \param mass robot mass [kg]
     */
    Eigen::Vector6d toState(double mass) const;
  };

  /** \brief Reference data. */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector2d pos = Eigen::Vector2d::Zero();

    //! CoM linear velocity [m/s]
    Eigen::Vector2d vel = Eigen::Vector2d::Zero();

    //! Angular momentum [kg m^2/s]
    Eigen::Vector2d angular_momentum = Eigen::Vector2d::Zero();

    /** \brief Get reference output size. */
    static constexpr int outputDim()
    {
      return 6;
    }

    /** \brief Get reference output of state-space model.
        \param mass robot mass [kg]
     */
    Eigen::Vector6d toOutput(double mass) const;
  };

  /** \brief Weight parameter. */
  struct WeightParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Linear momentum integral weight
    Eigen::Vector2d linear_momentum_integral;

    //! Linear momentum weight
    Eigen::Vector2d linear_momentum;

    //! Angular momentum weight
    Eigen::Vector2d angular_momentum;

    //! Force weight
    double force;

    /** \brief Constructor.
        \param _linear_momentum_integral linear momentum integral weight
        \param _linear_momentum linear momentum weight
        \param _angular_momentum angular momentum weight
        \param _force force weight
     */
    WeightParam(const Eigen::Vector2d & _linear_momentum_integral = Eigen::Vector2d::Constant(1.0),
                const Eigen::Vector2d & _linear_momentum = Eigen::Vector2d::Constant(0.0),
                const Eigen::Vector2d & _angular_momentum = Eigen::Vector2d::Constant(1.0),
                double _force = 1e-5)
    : linear_momentum_integral(_linear_momentum_integral), linear_momentum(_linear_momentum),
      angular_momentum(_angular_momentum), force(_force)
    {
    }

    /** \brief Get input weight vector.
        \param total_input_dim total input dimension
    */
    Eigen::VectorXd inputWeight(int total_input_dim) const;

    /** \brief Get output weight vector.
        \param seq_len sequence length
    */
    Eigen::VectorXd outputWeight(size_t seq_len) const;
  };

  /** \brief Motion data. */
  struct MotionData : MotionDataBase<Eigen::Vector2d, Eigen::Vector2d, Eigen::VectorXd>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Reference angular momentum [kg m^2/s]
    Eigen::Vector2d ref_angular_momentum;

    //! Planned angular momentum [kg m^2/s]
    Eigen::Vector2d planned_angular_momentum;

    /** \brief Dump data.
        \tparam StreamType stream type
     */
    template<class StreamType>
    void dump(StreamType & ofs) const
    {
      ofs << time << " " << ref_pos.transpose() << " " << ref_vel.transpose() << " " << ref_angular_momentum.transpose()
          << " " << planned_pos.transpose() << " " << planned_vel.transpose() << " " << planned_acc.transpose() << " "
          << planned_angular_momentum.transpose() << " " << planned_force.transpose() << std::endl;
    }
  };

  /** \brief State-space model. */
  class Model : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
        \param motion_param motion parameter
        \param output_dim output dimension
     */
    Model(double mass, const MotionParam & motion_param, int output_dim = 0);

  public:
    //! Motion parameter
    MotionParam motion_param_;
  };

  /** \brief State-space model for simulation. */
  class SimModel : public Model
  {
  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
        \param motion_param motion parameter
     */
    SimModel(double mass, const MotionParam & motion_param);
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param mass robot mass [kg]
      \param horizon_dt discretization timestep in horizon [sec]
      \param qp_solver_type QP solver type
  */
  LinearMpcXY(double mass,
              double horizon_dt,
              QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::QLD);

  /** \brief Plan one step.
      \param motion_param_func function of motion parameter
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param horizon_time_range start and end time of horizon ([sec], [sec])
      \param weight_param objective weight parameter
      \returns planned force sequence
  */
  Eigen::VectorXd planOnce(const std::function<MotionParam(double)> & motion_param_func,
                           const std::function<RefData(double)> & ref_data_func,
                           const InitialParam & initial_param,
                           const std::pair<double, double> & horizon_time_range,
                           const WeightParam & weight_param = WeightParam());

  /** \brief Plan with loop.
      \param motion_param_func function of motion parameter
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param motion_time_range start and end time of motion ([sec], [sec])
      \param horizon_duration horizon duration [sec]
      \param sim_dt discretization timestep for simulation [sec]
      \param weight_param objective weight parameter
      \returns planned force sequence
  */
  void planLoop(const std::function<MotionParam(double)> & motion_param_func,
                const std::function<RefData(double)> & ref_data_func,
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

  /** \brief Get motion data sequence. */
  inline const std::map<double, MotionData> & motionDataSeq() const
  {
    return motion_data_seq_;
  }

protected:
  /** \brief Process one step. */
  Eigen::VectorXd procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                           const StateDimVector & current_x,
                           const Eigen::VectorXd & ref_output_seq,
                           const WeightParam & weight_param);

public:
  //! Robot mass [kg]
  double mass_ = 0;

  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;

  //! Min/max z-component force [N]
  std::pair<double, double> force_range_;

protected:
  //! Motion data sequence
  std::map<double, MotionData> motion_data_seq_;
};
} // namespace CCC
