/* Author: Masaki Murooka */

#pragma once

#include <nmpc_ddp/DDPSolver.h>

namespace ForceColl
{
class Contact;
}

namespace CCC
{
/** \brief Differential dynamic programming (DDP) for centroidal model
           with single rigid-body dynamics (SRBD) approximation.

    Base link orientation is expressed in ZYX Euler angles.

    \todo Support other Euler angles
*/
class DdpSingleRigidBody
{
public:
  /** \brief Motion parameter. */
  struct MotionParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** \brief Contact list. */
    std::vector<std::shared_ptr<ForceColl::Contact>> contact_list;

    /** \brief Inertia matrix [kg m^2]. */
    Eigen::Matrix3d inertia_mat = Eigen::Matrix3d::Identity();
  };

  /** \brief Reference data. */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();

    //! Base link orientation [rad]
    Eigen::Vector3d ori = Eigen::Vector3d::Zero();
  };

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! CoM position weight in running cost
    Eigen::Vector3d running_pos;

    //! Base link orientation weight in running cost
    Eigen::Vector3d running_ori;

    //! Linear velocity weight in running cost
    Eigen::Vector3d running_linear_vel;

    //! Angular velocity weight in running cost
    Eigen::Vector3d running_angular_vel;

    //! Force weight in running cost
    double running_force;

    //! CoM position weight in terminal cost
    Eigen::Vector3d terminal_pos;

    //! Base link orientation weight in terminal cost
    Eigen::Vector3d terminal_ori;

    //! Linear velocity weight in terminal cost
    Eigen::Vector3d terminal_linear_vel;

    //! Angular velocity weight in terminal cost
    Eigen::Vector3d terminal_angular_vel;

    /** \brief Constructor.
        \param _running_pos CoM position weight in running cost
        \param _running_ori base link orientation weight in running cost
        \param _running_linear_vel linear velocity weight in running cost
        \param _running_angular_vel angular velocity weight in running cost
        \param _running_force force weight in running cost
        \param _terminal_pos CoM position weight in terminal cost
        \param _terminal_ori base link orientation weight in terminal cost
        \param _terminal_linear_vel linear velocity weight in terminal cost
        \param _terminal_angular_vel angular velocity weight in terminal cost
     */
    WeightParam(const Eigen::Vector3d & _running_pos = Eigen::Vector3d::Constant(1.0),
                const Eigen::Vector3d & _running_ori = Eigen::Vector3d::Constant(1.0),
                const Eigen::Vector3d & _running_linear_vel = Eigen::Vector3d::Constant(0.01),
                const Eigen::Vector3d & _running_angular_vel = Eigen::Vector3d::Constant(0.01),
                double _running_force = 1e-6,
                const Eigen::Vector3d & _terminal_pos = Eigen::Vector3d::Constant(1.0),
                const Eigen::Vector3d & _terminal_ori = Eigen::Vector3d::Constant(1.0),
                const Eigen::Vector3d & _terminal_linear_vel = Eigen::Vector3d::Constant(0.01),
                const Eigen::Vector3d & _terminal_angular_vel = Eigen::Vector3d::Constant(0.01))
    : running_pos(_running_pos), running_ori(_running_ori), running_linear_vel(_running_linear_vel),
      running_angular_vel(_running_angular_vel), running_force(_running_force), terminal_pos(_terminal_pos),
      terminal_ori(_terminal_ori), terminal_linear_vel(_terminal_linear_vel),
      terminal_angular_vel(_terminal_angular_vel)
    {
    }
  };

  /** \brief DDP problem of centroidal model with single rigid-body dynamics (SRBD) approximation.

      See #stateEq for the state equation of this problem.
   */
  class DdpProblem : public nmpc_ddp::DDPProblem<12, Eigen::Dynamic>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** \brief Constructor.
        \param horizon_dt discretization timestep in horizon [sec]
        \param mass robot mass [Kg]
        \param weight_param objective weight parameter
    */
    DdpProblem(double horizon_dt, double mass, const WeightParam & weight_param)
    : nmpc_ddp::DDPProblem<12, Eigen::Dynamic>(horizon_dt), mass_(mass), weight_param_(weight_param)
    {
    }

    /** \brief Set function of motion parameter.
        \param motion_param_func function of motion parameter
     */
    inline void setMotionParamFunc(const std::function<MotionParam(double)> & motion_param_func)
    {
      motion_param_func_ = motion_param_func;
    }

    /** \brief Set function of reference data.
        \param ref_data_func function of reference data
     */
    inline void setRefDataFunc(const std::function<RefData(double)> & ref_data_func)
    {
      ref_data_func_ = ref_data_func;
    }

    using DDPProblem::inputDim;

    /** \brief Gets the input dimension.
        \param t time
    */
    virtual int inputDim(double t) const override;

    /** \brief Calculate discrete state equation.
        \param t time [sec]
        \param x current state \f$\boldsymbol{x}_k\f$
        \param u current input \f$\boldsymbol{u}_k\f$
        \returns next state \f$\boldsymbol{x}_{k+1}\f$

        Dynamics is expressed by the following equation.
        \f{align*}{
        m \boldsymbol{\ddot{c}} &= \sum_i \lambda_i \boldsymbol{\rho}_i - m \boldsymbol{g} \\
        \boldsymbol{\dot{\omega}} + \boldsymbol{\omega} \times \boldsymbol{I} \boldsymbol{\omega} &= \sum_i
        (\boldsymbol{p}_i - \boldsymbol{c}) \times \lambda_i \boldsymbol{\rho}_i \f}

        This can be represented as a nonlinear system as follows.
        \f{align*}{
        \boldsymbol{\dot{x}} &=
        \begin{bmatrix}
          \boldsymbol{v} \\
          \boldsymbol{K}_{\mathit{Euler}}(\boldsymbol{\alpha}) \, \boldsymbol{\omega} \\
          \dfrac{1}{m} \sum_i \lambda_i \boldsymbol{\rho}_i - \boldsymbol{g} \\
          \boldsymbol{I}^{-1} \left( - \boldsymbol{\omega} \times \boldsymbol{I} \boldsymbol{\omega} +
          \sum_i (\boldsymbol{p}_i - \boldsymbol{c}) \times \lambda_i \boldsymbol{\rho}_i \right)
        \end{bmatrix}
        \f}

        State and control input are expressed by the following equations.
        \f{align*}{
        \boldsymbol{x} = \begin{bmatrix} \boldsymbol{c} \\ \boldsymbol{\alpha} \\ \boldsymbol{v} \\ \boldsymbol{\omega}
       \end{bmatrix}, \boldsymbol{u} = \begin{bmatrix} \vdots \\ \lambda_i \\ \vdots \end{bmatrix} \f}

        \f$m\f$ and \f$\boldsymbol{I}\f$ are the robot mass and inertia matrix, respectively.
        \f$\boldsymbol{c}\f$, \f$\boldsymbol{v}\f$, \f$\boldsymbol{\alpha}\f$, and \f$\boldsymbol{\omega}\f$ are
        CoM position, CoM velocity, base link orientation, and base link angular velocity, respectively.
        Base link orientation is expressed in Euler angles.
        \f$\boldsymbol{p}_i\f$, \f$\lambda_i\f$, and \f$\boldsymbol{\rho}_i\f$ are position, force scale, and ridge
        vector of i-th contact vertex ridge, respectively.

        Euler method is used to discretize the system.
        \f{align*}{
        \boldsymbol{x}_{k+1} = \boldsymbol{\dot{x}}_{k} \Delta t + \boldsymbol{x}_{k}
        \f}
     */
    virtual StateDimVector stateEq(double t, const StateDimVector & x, const InputDimVector & u) const override;

    /** \brief Calculate running cost.
        \param t time [sec]
        \param x current state \f$\boldsymbol{x}_k\f$
        \param u current input \f$\boldsymbol{u}_k\f$
        \returns running cost \f$L_k\f$
    */
    virtual double runningCost(double t, const StateDimVector & x, const InputDimVector & u) const override;

    /** \brief Calculate terminal cost.
        \param t time [sec]
        \param x current state \f$\boldsymbol{x}_k\f$
        \returns terminal cost \f$\phi_k\f$
    */
    virtual double terminalCost(double t, const StateDimVector & x) const override;

    /** \brief Calculate first-order derivatives of discrete state equation.
        \param t time [sec]
        \param x state
        \param u input
        \param state_eq_deriv_x first-order derivative of state equation w.r.t. state
        \param state_eq_deriv_u first-order derivative of state equation w.r.t. input

        The first-order derivative of the discrete state equation is expressed as follows.
        \f{align*}{
        \frac{\partial \boldsymbol{x}_{k+1}}{\partial \boldsymbol{x}_k} &=
        \begin{bmatrix}
          \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{E} & \boldsymbol{O} \\
          \boldsymbol{O} & \dfrac{\partial \boldsymbol{\dot{\alpha}}}{\partial \boldsymbol{\alpha}} & \boldsymbol{O} &
          \boldsymbol{K}_{\mathit{Euler}} \\
          \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} \\
          \boldsymbol{I}^{-1} (\sum_i \lambda_i \boldsymbol{\rho}_i) \times & \boldsymbol{O} & \boldsymbol{O} &
          \dfrac{\partial \boldsymbol{\dot{\omega}}}{\partial \boldsymbol{\omega}}
        \end{bmatrix} \Delta t + \boldsymbol{I} \\
        \frac{\partial \boldsymbol{x}_{k+1}}{\partial \boldsymbol{u}_k} &=
        \begin{bmatrix}
          \cdots & \boldsymbol{O} & \cdots \\
          \cdots & \boldsymbol{O} & \cdots \\
          \cdots & \dfrac{1}{m} \boldsymbol{\rho}_i & \cdots \\
          \cdots & \boldsymbol{I}^{-1} (\boldsymbol{p}_i - \boldsymbol{c}) \times \boldsymbol{\rho}_i & \cdots
        \end{bmatrix} \Delta t
        \f}

        \f$\boldsymbol{E}\f$ is the identity matrix.
        The formulas for \f$\dfrac{\partial \boldsymbol{\dot{\alpha}}}{\partial \boldsymbol{\alpha}}\f$ and
        \f$\dfrac{\partial \boldsymbol{\dot{\omega}}}{\partial \boldsymbol{\omega}}\f$ are complex,
        so they were derived using the symbolic mathematics library (SymPy).
    */
    virtual void calcStateEqDeriv(double t,
                                  const StateDimVector & x,
                                  const InputDimVector & u,
                                  Eigen::Ref<StateStateDimMatrix> state_eq_deriv_x,
                                  Eigen::Ref<StateInputDimMatrix> state_eq_deriv_u) const override;

    /** \brief Calculate first-order and second-order derivatives of discrete state equation.
        \param t time [sec]
        \param x state
        \param u input
        \param state_eq_deriv_x first-order derivative of state equation w.r.t. state
        \param state_eq_deriv_u first-order derivative of state equation w.r.t. input
        \param state_eq_deriv_xx second-order derivative of state equation w.r.t. state
        \param state_eq_deriv_uu second-order derivative of state equation w.r.t. input
        \param state_eq_deriv_xu second-order derivative of state equation w.r.t. state and input
    */
    inline virtual void calcStateEqDeriv(double, // t
                                         const StateDimVector &, // x
                                         const InputDimVector &, // u
                                         Eigen::Ref<StateStateDimMatrix>, // state_eq_deriv_x
                                         Eigen::Ref<StateInputDimMatrix>, // state_eq_deriv_u
                                         std::vector<StateStateDimMatrix> &, // state_eq_deriv_xx
                                         std::vector<InputInputDimMatrix> &, // state_eq_deriv_uu
                                         std::vector<StateInputDimMatrix> & // state_eq_deriv_xu
    ) const override
    {
      throw std::runtime_error("Second-order derivatives of state equation are not implemented.");
    }

    /** \brief Calculate first-order derivatives of running cost.
        \param t time [sec]
        \param x state
        \param u input
        \param running_cost_deriv_x first-order derivative of running cost w.r.t. state
        \param running_cost_deriv_u first-order derivative of running cost w.r.t. input
    */
    virtual void calcRunningCostDeriv(double t,
                                      const StateDimVector & x,
                                      const InputDimVector & u,
                                      Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                      Eigen::Ref<InputDimVector> running_cost_deriv_u) const override;

    /** \brief Calculate first-order and second-order derivatives of running cost.
        \param t time [sec]
        \param x state
        \param u input
        \param running_cost_deriv_x first-order derivative of running cost w.r.t. state
        \param running_cost_deriv_u first-order derivative of running cost w.r.t. input
        \param running_cost_deriv_xx second-order derivative of running cost w.r.t. state
        \param running_cost_deriv_uu second-order derivative of running cost w.r.t. input
        \param running_cost_deriv_xu second-order derivative of running cost w.r.t. state and input
    */
    virtual void calcRunningCostDeriv(double t,
                                      const StateDimVector & x,
                                      const InputDimVector & u,
                                      Eigen::Ref<StateDimVector> running_cost_deriv_x,
                                      Eigen::Ref<InputDimVector> running_cost_deriv_u,
                                      Eigen::Ref<StateStateDimMatrix> running_cost_deriv_xx,
                                      Eigen::Ref<InputInputDimMatrix> running_cost_deriv_uu,
                                      Eigen::Ref<StateInputDimMatrix> running_cost_deriv_xu) const override;

    /** \brief Calculate first-order derivatives of terminal cost.
        \param t time [sec]
        \param x state
        \param terminal_cost_deriv_x first-order derivative of terminal cost w.r.t. state
    */
    virtual void calcTerminalCostDeriv(double t,
                                       const StateDimVector & x,
                                       Eigen::Ref<StateDimVector> terminal_cost_deriv_x) const override;

    /** \brief Calculate first-order and second-order derivatives of terminal cost.
        \param t time [sec]
        \param x state
        \param terminal_cost_deriv_x first-order derivative of terminal cost w.r.t. state
        \param terminal_cost_deriv_xx second-order derivative of terminal cost w.r.t. state
    */
    virtual void calcTerminalCostDeriv(double t,
                                       const StateDimVector & x,
                                       Eigen::Ref<StateDimVector> terminal_cost_deriv_x,
                                       Eigen::Ref<StateStateDimMatrix> terminal_cost_deriv_xx) const override;

  public:
    //! Robot mass [Kg]
    double mass_ = 0;

  protected:
    //! Weight parameter
    WeightParam weight_param_;

    //! Function of motion parameter
    std::function<MotionParam(double)> motion_param_func_;

    //! Function of reference data
    std::function<RefData(double)> ref_data_func_;
  };

  /** \brief Initial parameter. */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();

    //! Base link orientation [rad]
    Eigen::Vector3d ori = Eigen::Vector3d::Zero();

    //! Linear velocity [m/s]
    Eigen::Vector3d linear_vel = Eigen::Vector3d::Zero();

    //! Angular velocity [rad/s]
    Eigen::Vector3d angular_vel = Eigen::Vector3d::Zero();

    /** \brief Initial guess of input sequence.

        Sequence length should be horizon_steps.
        If empty, initial guess of input sequence is initialized to all zeros.
    */
    std::vector<DdpProblem::InputDimVector> u_list = {};

    /** \brief Constructor. */
    InitialParam() {}

    /** \brief Constructor.
        \param state state of DDP problem
    */
    InitialParam(const DdpProblem::StateDimVector & state);

    /** \brief Get state of DDP problem. */
    DdpProblem::StateDimVector toState() const;
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param mass robot mass [kg]
      \param horizon_dt discretization timestep in horizon [sec]
      \param horizon_steps number of steps in horizon
      \param weight_param objective weight parameter
   */
  DdpSingleRigidBody(double mass,
                     double horizon_dt,
                     int horizon_steps,
                     const WeightParam & weight_param = WeightParam());

  /** \brief Plan one step.
      \param motion_param_func function of motion parameter
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned force scales
   */
  Eigen::VectorXd planOnce(const std::function<MotionParam(double)> & motion_param_func,
                           const std::function<RefData(double)> & ref_data_func,
                           const InitialParam & initial_param,
                           double current_time);

public:
  //! DDP problem
  std::shared_ptr<DdpProblem> ddp_problem_;

  //! DDP solver
  std::shared_ptr<nmpc_ddp::DDPSolver<12, Eigen::Dynamic>> ddp_solver_;

  //! Force scale limits (i.e., limits of \f$\lambda_i\f$ in the order of lower, upper)
  std::array<double, 2> force_scale_limits_ = {0.0, 1e6};
};
} // namespace CCC
