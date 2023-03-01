/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolverCollection.h>

#include <CCC/EigenTypes.h>
#include <CCC/VariantSequentialExtension.h>

namespace CCC
{
/** \brief Linear MPC of translational and rotational x,y-component motion considering predefined z-component motion.

    See the following for a detailed formulation.
      - H Audren, et al. Model preview control in multi-contact motion-application to a humanoid robot. IROS, 2014.
      - 長阪憲一郎, et al.
   接触拘束を考慮可能なマルチコンタクト対応スタビライザと一般化逆動力学による人型ロボットの全身制御.
   ロボティクスシンポジア予稿集, 2012.
 */
class LinearMpcXY
{
public:
  /** \brief State dimension. */
  static constexpr int state_dim_ = 6;

  /** \brief Type of state-space model. */
  using _StateSpaceModel = StateSpaceModel<state_dim_, Eigen::Dynamic, Eigen::Dynamic>;

  /** \brief Type of state vector. */
  using StateDimVector = _StateSpaceModel::StateDimVector;

public:
  /** \brief Motion parameter. */
  struct MotionParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM z position [m]
    double com_z;

    //! Total z force [N]
    double total_force_z;

    /** \brief List of contact vertex and force direction (i.e., friction pyramid ridge)

        Each column represents a pair of vertex and ridge. The first three rows represent vertex and the last three rows
       represent ridge.
     */
    Eigen::Matrix<double, 6, Eigen::Dynamic> vertex_ridge_list;
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

  /** \brief State-space model.

      Dynamics is expressed by the following equation.
      \f{align*}{
      \boldsymbol{\dot{P}} &= \sum_i \lambda_i \boldsymbol{\rho}_i - m \boldsymbol{g} \\
      \boldsymbol{\dot{L}} &= \sum_i (\boldsymbol{p}_i - \boldsymbol{c}) \times \lambda_i \boldsymbol{\rho}_i
      \f}
      \f$\boldsymbol{c}\f$, \f$\boldsymbol{P}\f$, and \f$\boldsymbol{L}\f$ are CoM, linear momentum, and angular
     momentum, respectively. \f$\boldsymbol{p}_i\f$, \f$\lambda_i\f$, and \f$\boldsymbol{\rho}_i\f$ are position, force
     scale, and ridge vector of i-th contact vertex ridge, respectively.

      This can be represented as a linear time-variant system as follows.
      \f{align*}{
      \boldsymbol{\dot{x}} &=
      \begin{bmatrix}
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 1 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & - \dfrac{f_z}{m} & 0 & 0 & 0 \\
        \dfrac{f_z}{m} & 0 & 0 & 0 & 0 & 0
      \end{bmatrix}
      \boldsymbol{x} +
      \begin{bmatrix}
        \cdots & 0 & \cdots \\
        \cdots & \rho_{i,x} & \cdots \\
        \cdots & 0 & \cdots \\
        \cdots & \rho_{i,y} & \cdots \\
        \cdots & - (p_{i,z} - c_z) \rho_{i,y} + p_{i,y} \rho_{i,z} & \cdots \\
        \cdots & (p_{i,z} - c_z) \rho_{i,x} - p_{i,x} \rho_{i,z} & \cdots
      \end{bmatrix} \boldsymbol{u}
      \f}
      We assume that vertical moton (i.e., CoM height \f$c_z\f$ and total vertical contact force \f$f_z\f$) is
     pre-defined.

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} m c_x \\ P_x \\ m c_y \\ P_y \\ L_x \\ L_y \end{bmatrix},
      \boldsymbol{u} = \begin{bmatrix} \vdots \\ \lambda_i \\ \vdots \end{bmatrix}
      \f}
   */
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
              QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any);

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

  //! Min/max ridge force [N]
  std::pair<double, double> force_range_;
};
} // namespace CCC
