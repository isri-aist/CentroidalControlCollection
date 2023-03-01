/* Author: Masaki Murooka */

#pragma once

#include <SpaceVecAlg/SpaceVecAlg>

#include <CCC/Constants.h>
#include <CCC/StateSpaceModel.h>

/** \brief State-space model of one-dimensional CoM-ZMP dynamics with ZMP input. */
class ComZmpSimModel1d : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param sim_dt simulation timestep [sec]
  */
  ComZmpSimModel1d(double com_height, double sim_dt)
  {
    double omega2 = CCC::constants::g / com_height;

    A_(0, 1) = 1;
    A_(1, 0) = omega2;

    B_(1, 0) = -1 * omega2;

    calcDiscMatrix(sim_dt);
  }

  /** \brief Update.
      \param input simulation input (ZMP)
  */
  void update(double input)
  {
    state_ = stateEqDisc(state_, (InputDimVector() << input).finished());
  }

public:
  //! Simulation state
  StateDimVector state_ = StateDimVector::Zero();
};

/** \brief State-space model of vertical motion dynamics with force input. */
class VerticalSimModel : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param mass robot mass [Kg]
      \param sim_dt simulation timestep [sec]
  */
  VerticalSimModel(double mass, double sim_dt)
  {
    A_(0, 1) = 1;

    B_(1, 0) = 1 / mass;

    E_(1) = -1 * CCC::constants::g;

    calcDiscMatrix(sim_dt);
  }

  /** \brief Update.
      \param input simulation input (force)
  */
  void update(double input)
  {
    state_ = stateEqDisc(state_, (InputDimVector() << input).finished());
  }

public:
  //! Simulation state
  StateDimVector state_ = StateDimVector::Zero();
};

/** \brief Simulation of two-dimensional CoM-ZMP dynamics with ZMP input. */
class ComZmpSim2d
{
public:
  /** \brief Simulation state. */
  struct State
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! X position and velocity
    Eigen::Vector2d x = Eigen::Vector2d::Zero();

    //! Y position and velocity
    Eigen::Vector2d y = Eigen::Vector2d::Zero();

    /** \brief Get position. */
    inline Eigen::Vector2d pos() const
    {
      return Eigen::Vector2d(x[0], y[0]);
    }

    /** \brief Get velocity. */
    inline Eigen::Vector2d vel() const
    {
      return Eigen::Vector2d(x[1], y[1]);
    }
  };

  /** \brief Simulation input (ZMP). */
  using Input = Eigen::Vector2d;

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param sim_dt simulation timestep [sec]
  */
  ComZmpSim2d(double com_height, double sim_dt) : model_(std::make_shared<ComZmpSimModel1d>(com_height, sim_dt)) {}

  /** \brief Update.
      \param input simulation input (ZMP)
  */
  void update(const Input & input)
  {
    state_.x = model_->stateEqDisc(state_.x, (ComZmpSimModel1d::InputDimVector() << input.x()).finished());
    state_.y = model_->stateEqDisc(state_.y, (ComZmpSimModel1d::InputDimVector() << input.y()).finished());
  }

  /** \brief Add disturbance
      \param impulse_per_mass impulse per mass [m/s]
  */
  void addDisturb(const Eigen::Vector2d & impulse_per_mass)
  {
    state_.x[1] += impulse_per_mass.x();
    state_.y[1] += impulse_per_mass.x();
  }

public:
  //! One-dimensional CoM-ZMP model
  std::shared_ptr<ComZmpSimModel1d> model_;

  //! Simulation state
  State state_;
};

/** \brief Simulation of three-dimensional CoM-ZMP dynamics with ZMP and force input. */
class ComZmpSim3d
{
public:
  /** \brief Simulation state. */
  struct State
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! X position and velocity
    Eigen::Vector2d x = Eigen::Vector2d::Zero();

    //! Y position and velocity
    Eigen::Vector2d y = Eigen::Vector2d::Zero();

    //! Z position and velocity
    Eigen::Vector2d z = Eigen::Vector2d::Zero();

    /** \brief Get position. */
    inline Eigen::Vector3d pos() const
    {
      return Eigen::Vector3d(x[0], y[0], z[0]);
    }

    /** \brief Get velocity. */
    inline Eigen::Vector3d vel() const
    {
      return Eigen::Vector3d(x[1], y[1], z[1]);
    }
  };

  /** \brief Simulation input. */
  struct Input
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! ZMP [m]
    Eigen::Vector2d zmp = Eigen::Vector2d::Zero();

    //! Z force [N]
    double force_z = 0;
  };

public:
  /** \brief Constructor.
      \param mass robot mass [Kg]
      \param sim_dt simulation timestep [sec]
  */
  ComZmpSim3d(double mass, double sim_dt) : sim_dt_(sim_dt), z_model_(std::make_shared<VerticalSimModel>(mass, sim_dt))
  {
  }

  /** \brief Update.
      \param input simulation input (ZMP)
  */
  void update(const Input & input)
  {
    xy_model_ = std::make_shared<ComZmpSimModel1d>(state_.z[0], sim_dt_);

    state_.x = xy_model_->stateEqDisc(state_.x, (ComZmpSimModel1d::InputDimVector() << input.zmp.x()).finished());
    state_.y = xy_model_->stateEqDisc(state_.y, (ComZmpSimModel1d::InputDimVector() << input.zmp.y()).finished());
    state_.z = z_model_->stateEqDisc(state_.z, (VerticalSimModel::InputDimVector() << input.force_z).finished());
  }

  /** \brief Add disturbance
      \param impulse_per_mass impulse per mass [m/s]
  */
  void addDisturb(const Eigen::Vector2d & impulse_per_mass)
  {
    state_.x[1] += impulse_per_mass.x();
    state_.y[1] += impulse_per_mass.x();
  }

public:
  //! Simulation timestep [sec]
  double sim_dt_ = 0;

  //! One-dimensional CoM-ZMP model
  std::shared_ptr<ComZmpSimModel1d> xy_model_;

  //! Vertical motion model
  std::shared_ptr<VerticalSimModel> z_model_;

  //! Simulation state
  State state_;
};

/** \brief Simulation of centroidal dynamics.

    The following approximations are made for rotational motion.
    - The diagonal elements of the moment of inertia are constant, and the off-diagonal elements are zero.
    - Angular velocity and the derivatives of Euler angle are considered the same.
    - Ignore the nonlinear term of angular momentum.
 */
class CentroidalSim : public CCC::StateSpaceModel<18, 6, 0>
{
public:
  /** \brief Simulation state. */
  struct State
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m] and base link orientation [rad]
    sva::MotionVecd pos = sva::MotionVecd::Zero();

    //! CoM velocity [m/s] and base link angular velocity [rad/s]
    sva::MotionVecd vel = sva::MotionVecd::Zero();

    //! Linear momentum [kg m/s] and angular momentum [kg m^2/s]
    sva::ForceVecd momentum = sva::ForceVecd::Zero();

    /** \brief Constructor. */
    State() {}

    /** \brief Constructor.
        \param state_vec state of state-space model
     */
    State(const StateDimVector & state_vec)
    {
      pos = sva::MotionVecd(state_vec.segment<3>(3), state_vec.segment<3>(0));
      vel = sva::MotionVecd(state_vec.segment<3>(9), state_vec.segment<3>(6));
      momentum = sva::ForceVecd(state_vec.segment<3>(15), state_vec.segment<3>(12));
    }

    /** \brief Get state of state-space model. */
    inline StateDimVector toState() const
    {
      StateDimVector state_vec;
      state_vec << pos.linear(), pos.angular(), vel.linear(), vel.angular(), momentum.force(), momentum.moment();
      return state_vec;
    }
  };

  /** \brief Simulation input.

      Simulation input is wrench only.
   */
  struct Input : public sva::ForceVecd
  {
    /** \brief Constructor. */
    Input() {}

    /** \brief Constructor.
        \param wrench input wrench
     */
    Input(const sva::ForceVecd & wrench) : sva::ForceVecd(wrench) {}

    /** \brief Get input of state-space model. */
    inline InputDimVector toInput() const
    {
      InputDimVector input_vec;
      input_vec << force(), moment();
      return input_vec;
    }
  };

public:
  /** \brief Constructor.
      \param mass robot mass [Kg]
      \param moment_of_inertia moment of inertia of robot [kg m^2]
      \param sim_dt simulation timestep [sec]
  */
  CentroidalSim(double mass, const Eigen::Vector3d & moment_of_inertia, double sim_dt) : sim_dt_(sim_dt)
  {
    A_.block<6, 6>(0, 6).diagonal().setConstant(1.0);

    B_.block<3, 3>(6, 0).diagonal().setConstant(1.0 / mass);
    B_.block<3, 3>(9, 3).diagonal() = moment_of_inertia.cwiseInverse();
    B_.block<6, 6>(12, 0).diagonal().setConstant(1.0);

    E_.segment<3>(6).z() = -1 * CCC::constants::g;
    E_.segment<3>(12).z() = -1 * mass * CCC::constants::g;

    calcDiscMatrix(sim_dt_);
  }

  /** \brief Update.
      \param input simulation input

      Moment is represented around CoM.
  */
  void update(const Input & input)
  {
    state_ = State(stateEqDisc(state_.toState(), input.toInput()));
  }

  /** \brief Add disturbance
      \param impulse_per_mass linear and angular impulse per mass [m/s], [m^2/s]
  */
  void addDisturb(const sva::ForceVecd & impulse_per_mass)
  {
    state_.vel.linear() += impulse_per_mass.force();
    state_.vel.angular() += impulse_per_mass.moment();
  }

public:
  //! Simulation timestep [sec]
  double sim_dt_ = 0;

  //! Simulation state
  State state_;
};
