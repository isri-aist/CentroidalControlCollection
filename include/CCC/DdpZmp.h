/* Author: Masaki Murooka */

#pragma once

#include <nmpc_ddp/DDPSolver.h>

namespace CCC
{
/** \brief Differential dynamic programming (DDP) for CoM-ZMP model.

    See the following for a detailed formulation.
      - S Feng, et al. Optimization‐based full body control for the darpa robotics challenge. Journal of field robotics,
   2015.
 */
class DdpZmp
{
public:
  /** \brief Reference data. */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! ZMP [m]
    Eigen::Vector3d zmp = Eigen::Vector3d::Zero();

    //! CoM z position [m]
    double com_z = 0;
  };

  /** \brief Planned data. */
  struct PlannedData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! ZMP [m]
    Eigen::Vector2d zmp = Eigen::Vector2d::Zero();

    //! Z force [N]
    double force_z = 0;
  };

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! CoM z position weight in running cost
    double running_com_pos_z;

    //! ZMP weight in running cost
    double running_zmp;

    //! Z force weight in running cost
    double running_force_z;

    //! CoM x and y position weight in terminal cost
    double terminal_com_pos_xy;

    //! CoM z position weight in terminal cost
    double terminal_com_pos_z;

    //! CoM velocity weight in terminal cost
    double terminal_com_vel;

    /** \brief Constructor.
        \param _running_com_pos_z CoM z weight in running cost
        \param _running_zmp ZMP weight in running cost
        \param _running_force_z z force weight in running cost
        \param _terminal_com_pos_xy CoM x and y weight in terminal cost
        \param _terminal_com_pos_z CoM z weight in terminal cost
        \param _terminal_com_vel CoM velocity weight in terminal cost
     */
    WeightParam(double _running_com_pos_z = 1e2,
                double _running_zmp = 1e-1,
                double _running_force_z = 1e-4,
                double _terminal_com_pos_xy = 1.0,
                double _terminal_com_pos_z = 1e2,
                double _terminal_com_vel = 1.0)
    : running_com_pos_z(_running_com_pos_z), running_zmp(_running_zmp), running_force_z(_running_force_z),
      terminal_com_pos_xy(_terminal_com_pos_xy), terminal_com_pos_z(_terminal_com_pos_z),
      terminal_com_vel(_terminal_com_vel)
    {
    }
  };

  /** \brief DDP problem of CoM-ZMP model.

      Dynamics is expressed by the following equation.
      \f{align*}{
      \begin{bmatrix} \ddot{c}_x \\ \ddot{c}_y \\ \ddot{c}_z \end{bmatrix} =
      \begin{bmatrix}
      \dfrac{(c_x - z_x) f_z}{m (c_z - z_z)} \\ \dfrac{(c_y - z_y) f_z}{m (c_z - z_z)} \\ \dfrac{f_z}{m} - g
      \end{bmatrix}
      \f}
      \f$\boldsymbol{c}\f$ and \f$\boldsymbol{z}\f$ are CoM and ZMP.

      State and control input are expressed by the following equations.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} c_x \\ \dot{c}_x \\ c_y \\ \dot{c}_y \\ c_z \\ \dot{c}_z \end{bmatrix},
      \boldsymbol{u} = \begin{bmatrix} z_x \\ z_y \\ f_z \end{bmatrix}
      \f}
   */
  class DdpProblem : public nmpc_ddp::DDPProblem<6, 3>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** \brief Constructor.
        \param horizon_dt discretization timestep in horizon [sec]
        \param mass robot mass [Kg]
        \param weight_param objective weight parameter
    */
    DdpProblem(double horizon_dt, double mass, const WeightParam & weight_param)
    : nmpc_ddp::DDPProblem<6, 3>(horizon_dt), mass_(mass), weight_param_(weight_param)
    {
    }

    /** \brief Set function of reference data.
        \param ref_data_func function of reference data
     */
    void setRefDataFunc(const std::function<RefData(double)> & ref_data_func)
    {
      ref_data_func_ = ref_data_func;
    }

    /** \brief Calculate discrete state equation.
        \param t time [sec]
        \param x current state (x[k])
        \param u current input (u[k])
        \returns next state (x[k+1])
     */
    virtual StateDimVector stateEq(double t, const StateDimVector & x, const InputDimVector & u) const override;

    /** \brief Calculate running cost.
        \param t time [sec]
        \param x current state (x[k])
        \param u current input (u[k])
        \returns running cost (L[k])
    */
    virtual double runningCost(double t, const StateDimVector & x, const InputDimVector & u) const override;

    /** \brief Calculate terminal cost.
        \param t time [sec]
        \param x current state (x[k])
        \returns terminal cost (phi[k])
    */
    virtual double terminalCost(double t, const StateDimVector & x) const override;

    /** \brief Calculate first-order derivatives of discrete state equation.
        \param t time [sec]
        \param x state
        \param u input
        \param state_eq_deriv_x first-order derivative of state equation w.r.t. state
        \param state_eq_deriv_u first-order derivative of state equation w.r.t. input
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

  protected:
    //! Robot mass [Kg]
    double mass_;

    //! Weight parameter
    WeightParam weight_param_;

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

    /** \brief Initial guess of input sequence.

        Sequence length should be horizon_steps.
        If empty, initial guess of input sequence is initialized to all zeros.
    */
    std::vector<DdpProblem::InputDimVector> u_list = {};

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
  DdpZmp(double mass, double horizon_dt, int horizon_steps, const WeightParam & weight_param = WeightParam())
  : ddp_problem_(std::make_shared<DdpProblem>(horizon_dt, mass, weight_param)),
    ddp_solver_(std::make_shared<nmpc_ddp::DDPSolver<6, 3>>(ddp_problem_))
  {
    ddp_solver_->config().horizon_steps = horizon_steps;
  }

  /** \brief Plan one step.
      \param ref_data_func function of reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned data
   */
  PlannedData planOnce(const std::function<RefData(double)> & ref_data_func,
                       const InitialParam & initial_param,
                       double current_time);

public:
  //! DDP problem
  std::shared_ptr<DdpProblem> ddp_problem_;

  //! DDP solver
  std::shared_ptr<nmpc_ddp::DDPSolver<6, 3>> ddp_solver_;
};
} // namespace CCC
