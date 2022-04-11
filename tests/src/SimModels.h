/* Author: Masaki Murooka */

#pragma once

#include <memory>

#include <CCC/Constants.h>
#include <CCC/StateSpaceModel.h>

/** \brief State-space model of one-dimensional CoM-ZMP dynamics with ZMP input. */
class ComZmpSimModel1d : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
  */
  ComZmpSimModel1d(double com_height)
  {
    double omega2 = CCC::constants::g / com_height;

    A_(0, 1) = 1;
    A_(1, 0) = omega2;

    B_(1, 0) = -1 * omega2;
  }
};

/** \brief State-space model of vertical motion dynamics with force input. */
class VerticalSimModel : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param mass robot mass [Kg]
  */
  VerticalSimModel(double mass)
  {
    A_(0, 1) = 1;

    B_(1, 0) = 1 / mass;

    E_(1) = -1 * CCC::constants::g;
  }
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
  ComZmpSim2d(double com_height, double sim_dt) : model_(std::make_shared<ComZmpSimModel1d>(com_height))
  {
    model_->calcDiscMatrix(sim_dt);
  }

  /** \brief Update.
      \param input simulation input (ZMP)
  */
  void update(const Input & input)
  {
    ComZmpSimModel1d::InputDimVector input_1d;
    input_1d << input.x();
    state_.x = model_->stateEqDisc(state_.x, input_1d);
    input_1d << input.y();
    state_.y = model_->stateEqDisc(state_.y, input_1d);
  }

public:
  //! One-dimensional CoM-ZMP model
  std::shared_ptr<ComZmpSimModel1d> model_;

  //! Simulation state
  State state_;
};
