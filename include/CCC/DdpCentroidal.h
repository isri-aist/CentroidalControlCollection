/* Author: Masaki Murooka */

#pragma once

#include <nmpc_ddp/DDPSolver.h>

#include <CCC/EigenTypes.h>
#include <CCC/MotionData.h>

#include <map>

namespace CCC
{
/** \brief Differential dynamic programming (DDP) for centroidal model. */
class DdpCentroidal
{
public:
  /** \brief Motion parameter. */
  struct MotionParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** \brief List of contact vertex and force direction (i.e., friction pyramid ridge)

        Each column represents a pair of vertex and ridge. The first three rows represent vertex and the last three rows
       represent ridge.
     */
    Eigen::Matrix<double, 6, Eigen::Dynamic> vertex_ridge_list;

    /** \brief Calculate total wrench.
        \param force_scales force scales
        \param moment_origin moment origin
        \returns wrench (in order of force, moment)
    */
    Eigen::Vector6d calcTotalWrench(const Eigen::VectorXd & force_scales,
                                    const Eigen::Vector3d & moment_origin = Eigen::Vector3d::Zero()) const;
  };

  /** \brief Reference data. */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  };

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! CoM position weight in running cost
    Eigen::Vector3d running_pos;

    //! Linear momentum weight in running cost
    Eigen::Vector3d running_linear_momentum;

    //! Angular momentum weight in running cost
    Eigen::Vector3d running_angular_momentum;

    //! Force weight in running cost
    double running_force;

    //! CoM position weight in terminal cost
    Eigen::Vector3d terminal_pos;

    //! Linear momentum weight in terminal cost
    Eigen::Vector3d terminal_linear_momentum;

    //! Angular momentum weight in terminal cost
    Eigen::Vector3d terminal_angular_momentum;

    /** \brief Constructor.
        \param _running_pos CoM position weight in running cost
        \param _running_linear_momentum linear momentum weight in running cost
        \param _running_angular_momentum angular momentum weight in running cost
        \param _running_force force weight in running cost
        \param _terminal_pos CoM position weight in terminal cost
        \param _terminal_linear_momentum linear momentum weight in terminal cost
        \param _terminal_angular_momentum angular momentum weight in terminal cost
     */
    WeightParam(const Eigen::Vector3d & _running_pos = Eigen::Vector3d::Constant(1.0),
                const Eigen::Vector3d & _running_linear_momentum = Eigen::Vector3d::Constant(0.0),
                const Eigen::Vector3d & _running_angular_momentum = Eigen::Vector3d::Constant(1.0),
                double _running_force = 1e-6,
                const Eigen::Vector3d & _terminal_pos = Eigen::Vector3d::Constant(1.0),
                const Eigen::Vector3d & _terminal_linear_momentum = Eigen::Vector3d::Constant(0.0),
                const Eigen::Vector3d & _terminal_angular_momentum = Eigen::Vector3d::Constant(1.0))
    : running_pos(_running_pos), running_linear_momentum(_running_linear_momentum),
      running_angular_momentum(_running_angular_momentum), running_force(_running_force), terminal_pos(_terminal_pos),
      terminal_linear_momentum(_terminal_linear_momentum), terminal_angular_momentum(_terminal_angular_momentum)
    {
    }
  };

  /** \brief Motion data. */
  struct MotionData : MotionDataBase<Eigen::Vector3d, Eigen::Vector3d, Eigen::VectorXd>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Reference angular momentum [kg m^2/s]
    Eigen::Vector3d ref_angular_momentum;

    //! Planned angular momentum [kg m^2/s]
    Eigen::Vector3d planned_angular_momentum;

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

  /** \brief DDP problem of centroidal model.

      See #stateEq for the state equation of this problem.
   */
  class DdpProblem : public nmpc_ddp::DDPProblem<9, Eigen::Dynamic>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** \brief Constructor.
        \param horizon_dt discretization timestep in horizon [sec]
        \param mass robot mass [Kg]
        \param weight_param objective weight parameter
    */
    DdpProblem(double horizon_dt, double mass, const WeightParam & weight_param)
    : nmpc_ddp::DDPProblem<9, Eigen::Dynamic>(horizon_dt), mass_(mass), weight_param_(weight_param)
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
    virtual inline int inputDim(double t) const override
    {
      const MotionParam & motion_param = motion_param_func_(t);
      return motion_param.vertex_ridge_list.cols();
    }

    /** \brief Calculate discrete state equation.
        \param t time [sec]
        \param x current state \f$\boldsymbol{x}_k\f$
        \param u current input \f$\boldsymbol{u}_k\f$
        \returns next state \f$\boldsymbol{x}_{k+1}\f$

        Dynamics is expressed by the following equation.
        \f{align*}{
        \boldsymbol{\dot{P}} &= \sum_i \lambda_i \boldsymbol{\rho}_i - m \boldsymbol{g} \\
        \boldsymbol{\dot{L}} &= \sum_i (\boldsymbol{p}_i - \boldsymbol{c}) \times \lambda_i \boldsymbol{\rho}_i
        \f}
        \f$\boldsymbol{c}\f$, \f$\boldsymbol{P}\f$, and \f$\boldsymbol{L}\f$ are CoM, linear momentum, and angular
       momentum, respectively. \f$\boldsymbol{p}_i\f$, \f$\lambda_i\f$, and \f$\boldsymbol{\rho}_i\f$ are position,
       force scale, and ridge vector of i-th contact vertex ridge, respectively.

        This can be represented as a nonlinear system as follows.
        \f{align*}{
        \boldsymbol{\dot{x}} &=
        \begin{bmatrix}
          \dfrac{1}{m} \boldsymbol{P} \\
          \sum_i \lambda_i \boldsymbol{\rho}_i - m \boldsymbol{g} \\
          \sum_i (\boldsymbol{p}_i - \boldsymbol{c}) \times \lambda_i \boldsymbol{\rho}_i
        \end{bmatrix}
        \f}

        State and control input are expressed by the following equations.
        \f{align*}{
        \boldsymbol{x} = \begin{bmatrix} \boldsymbol{c} \\ \boldsymbol{P} \\ \boldsymbol{L} \end{bmatrix},
        \boldsymbol{u} = \begin{bmatrix} \vdots \\ \lambda_i \\ \vdots \end{bmatrix}
        \f}

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
          \boldsymbol{O} & \dfrac{1}{m} \boldsymbol{I} & \boldsymbol{O} \\
          \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} \\
          (\sum_i \lambda_i \boldsymbol{\rho}_i) \times & \boldsymbol{O} & \boldsymbol{O}
        \end{bmatrix} \Delta t + \boldsymbol{I} \\
        \frac{\partial \boldsymbol{x}_{k+1}}{\partial \boldsymbol{u}_k} &=
        \begin{bmatrix}
          \cdots & \boldsymbol{O} & \cdots \\
          \cdots & \boldsymbol{\rho}_i & \cdots \\
          \cdots & (\boldsymbol{p}_i - \boldsymbol{c}) \times \boldsymbol{\rho}_i & \cdots
        \end{bmatrix} \Delta t
        \f}
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
    inline virtual void calcStateEqDeriv(double t,
                                         const StateDimVector & x,
                                         const InputDimVector & u,
                                         Eigen::Ref<StateStateDimMatrix> state_eq_deriv_x,
                                         Eigen::Ref<StateInputDimMatrix> state_eq_deriv_u,
                                         std::vector<StateStateDimMatrix> & state_eq_deriv_xx,
                                         std::vector<InputInputDimMatrix> & state_eq_deriv_uu,
                                         std::vector<StateInputDimMatrix> & state_eq_deriv_xu) const override
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

    //! CoM velocity [m/s]
    Eigen::Vector3d vel = Eigen::Vector3d::Zero();

    //! Angular momentum [kg m^2/s]
    Eigen::Vector3d angular_momentum = Eigen::Vector3d::Zero();

    /** \brief Initial guess of input sequence.

        Sequence length should be horizon_steps.
        If empty, initial guess of input sequence is initialized to all zeros.
    */
    std::vector<DdpProblem::InputDimVector> u_list = {};

    /** \brief Constructor. */
    InitialParam() {}

    /** \brief Constructor.
        \param state state of DDP problem
        \param mass robot mass [kg]
    */
    InitialParam(const DdpProblem::StateDimVector & state, double mass);

    /** \brief Get state of DDP problem.
        \param mass robot mass [kg]
     */
    DdpProblem::StateDimVector toState(double mass) const;
  };

  /** \brief State-space model for simulation. */
  class SimModel
  {
  public:
    /** \brief Type of state vector. */
    using StateDimVector = DdpProblem::StateDimVector;

    /** \brief Type of input vector. */
    using InputDimVector = DdpProblem::InputDimVector;

    /** \brief Type of output vector. */
    using OutputDimVector = Eigen::Matrix<double, 12, 1>;

  public:
    /** \brief Constructor.
        \param mass robot mass [kg]
        \param motion_param_func function of motion parameter
        \param dt discretization timestep [sec]
     */
    SimModel(double mass, const std::function<MotionParam(double)> & motion_param_func, double dt)
    : mass_(mass), motion_param_func_(motion_param_func), dt_(dt)
    {
    }

    /** \brief Calculate next state and observation output.
        \param x current state
        \param u current input
        \param next_x next state
        \param y current output
    */
    void procOnce(double t,
                  const StateDimVector & x,
                  const InputDimVector & u,
                  Eigen::Ref<StateDimVector> next_x,
                  Eigen::Ref<OutputDimVector> y) const;

  public:
    //! Robot mass [Kg]
    double mass_ = 0;

    //! Function of motion parameter
    std::function<MotionParam(double)> motion_param_func_;

    //! Discretization timestep [sec]
    double dt_ = 0;
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param mass robot mass [kg]
      \param horizon_dt discretization timestep in horizon [sec]
      \param horizon_steps number of steps in horizon
      \param weight_param objective weight parameter
   */
  DdpCentroidal(double mass, double horizon_dt, int horizon_steps, const WeightParam & weight_param = WeightParam());

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

  /** \brief Plan with loop.
      \param motion_param_func function of motion parameter
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param motion_time_range start and end time of motion ([sec], [sec])
      \param sim_dt discretization timestep for simulation [sec]
      \param ddp_max_iter DDP max iteration for the second and subsequent loop iterations (for the first loop iteration,
     the currently set DDP max iteration is used)
  */
  void planLoop(const std::function<MotionParam(double)> & motion_param_func,
                const std::function<RefData(double)> & ref_data_func,
                const InitialParam & initial_param,
                const std::pair<double, double> & motion_time_range,
                double sim_dt,
                int ddp_max_iter = 1);

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

public:
  //! DDP problem
  std::shared_ptr<DdpProblem> ddp_problem_;

  //! DDP solver
  std::shared_ptr<nmpc_ddp::DDPSolver<9, Eigen::Dynamic>> ddp_solver_;

  //! Force scale limits (i.e., limits of \f$\lambda_i\f$ in the order of lower, upper)
  std::array<double, 2> force_scale_limits_ = {0.0, 1e6};

protected:
  //! Motion data sequence
  std::map<double, MotionData> motion_data_seq_;
};
} // namespace CCC
