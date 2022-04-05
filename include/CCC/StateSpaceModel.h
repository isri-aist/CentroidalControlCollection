/* Author: Masaki Murooka */

#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace CCC
{
/** \brief State-space model.

    dx = A x + B u + E
    y = C x + D u + F

    \tparam StateDim state dimension
    \tparam InputDim input dimension
    \tparam OutputDim output dimension
*/
template<int StateDim, int InputDim, int OutputDim>
class StateSpaceModel
{
public:
  /** \brief Type of state vector. */
  using StateDimVector = Eigen::Matrix<double, StateDim, 1>;

  /** \brief Type of input vector. */
  using InputDimVector = Eigen::Matrix<double, InputDim, 1>;

  /** \brief Type of output vector. */
  using OutputDimVector = Eigen::Matrix<double, OutputDim, 1>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param state_dim state dimension
      \param input_dim input dimension
      \param output_dim output dimension
      \note dimensions in parameter can be omitted if a fixed value is given in the template value.
  */
  StateSpaceModel(int state_dim = StateDim, int input_dim = InputDim, int output_dim = OutputDim)
  : state_dim_(state_dim), input_dim_(input_dim), output_dim_(output_dim)
  {
    // Check dimension is positive
    if(state_dim_ <= 0)
    {
      throw std::runtime_error("state_dim must be positive: " + std::to_string(state_dim_) + " <= 0");
    }
    if(input_dim_ < 0)
    {
      throw std::runtime_error("input_dim must be non-negative: " + std::to_string(input_dim_) + " < 0");
    }
    if(output_dim_ < 0)
    {
      throw std::runtime_error("output_dim must be non-negative: " + std::to_string(output_dim_) + " < 0");
    }

    // Check dimension consistency
    if constexpr(StateDim != Eigen::Dynamic)
    {
      if(state_dim_ != StateDim)
      {
        throw std::runtime_error("state_dim is inconsistent with template parameter: " + std::to_string(state_dim_)
                                 + " != " + std::to_string(StateDim));
      }
    }
    if constexpr(InputDim != Eigen::Dynamic)
    {
      if(input_dim_ != InputDim)
      {
        throw std::runtime_error("input_dim is inconsistent with template parameter: " + std::to_string(input_dim_)
                                 + " != " + std::to_string(InputDim));
      }
    }
    if constexpr(OutputDim != Eigen::Dynamic)
    {
      if(output_dim_ != OutputDim)
      {
        throw std::runtime_error("output_dim is inconsistent with template parameter: " + std::to_string(output_dim_)
                                 + " != " + std::to_string(OutputDim));
      }
    }

    // Initialize with zero matrix if size is fixed
    A_.setZero(state_dim_, state_dim_);
    B_.setZero(state_dim_, input_dim_);
    C_.setZero(output_dim_, state_dim_);
    D_.setZero(output_dim_, input_dim_);
    E_.setZero(state_dim_);
    F_.setZero(output_dim_);
  }

  /** \brief Destructor.
      \note Need to make class polymorphic.
      See https://stackoverflow.com/a/15114118
  */
  virtual ~StateSpaceModel() = default;

  /** \brief Gets the state dimension. */
  int stateDim() const
  {
    return state_dim_;
  }

  /** \brief Gets the input dimension. */
  int inputDim() const
  {
    return input_dim_;
  }

  /** \brief Gets the output dimension. */
  int outputDim() const
  {
    return output_dim_;
  }

  /** \brief Calculate the continuous state equation.
      \param x state
      \param u input
      \returns time derivative of state (\dot{x})
  */
  StateDimVector stateEq(const StateDimVector & x, const InputDimVector & u) const
  {
    return A_ * x + B_ * u + E_;
  }

  /** \brief Calculate the discrete state equation.
      \param x current state (x[k])
      \param u current input (u[k])
      \returns next state (x[k+1])
  */
  StateDimVector stateEqDisc(const StateDimVector & x, const InputDimVector & u) const
  {
    return Ad_ * x + Bd_ * u + Ed_;
  }

  /** \brief Calculate the observation equation.
      \param x state
      \returns observation
  */
  OutputDimVector observEq(const StateDimVector & x) const
  {
    return C_ * x + F_;
  }

  /** \brief Calculate the observation equation.
      \param x state
      \param u input
      \returns observation
  */
  OutputDimVector observEq(const StateDimVector & x, const InputDimVector & u) const
  {
    return C_ * x + D_ * u + F_;
  }

  /** \brief Calculate the discrete system matrices.
      \param dt discretization timestep [sec]
  */
  void calcDiscMatrix(double dt)
  {
    dt_ = dt;

    // Zero-order hold discretization
    // See https://en.wikipedia.org/wiki/Discretization
    if constexpr(StateDim == Eigen::Dynamic || InputDim == Eigen::Dynamic)
    {
      // If dimension is dynamic
      if(E_.norm() == 0)
      {
        Eigen::MatrixXd ABZero(state_dim_ + input_dim_, state_dim_ + input_dim_);
        ABZero << dt_ * A_, dt_ * B_, Eigen::MatrixXd::Zero(input_dim_, state_dim_ + input_dim_);
        Eigen::MatrixXd ABZeroExp = ABZero.exp();
        Ad_ = ABZeroExp.template block(0, 0, state_dim_, state_dim_);
        Bd_ = ABZeroExp.template block(0, state_dim_, state_dim_, input_dim_);
        Ed_.setZero(state_dim_);
      }
      else
      {
        // There is no proof that this is correct.
        Eigen::MatrixXd ABEZero(state_dim_ + input_dim_ + 1, state_dim_ + input_dim_ + 1);
        ABEZero << dt_ * A_, dt_ * B_, dt_ * E_, Eigen::MatrixXd::Zero(input_dim_ + 1, state_dim_ + input_dim_ + 1);
        Eigen::MatrixXd ABEZeroExp = ABEZero.exp();
        Ad_ = ABEZeroExp.template block(0, 0, state_dim_, state_dim_);
        Bd_ = ABEZeroExp.template block(0, state_dim_, state_dim_, input_dim_);
        Ed_ = ABEZeroExp.template block(0, state_dim_ + input_dim_, state_dim_, 1);
      }
    }
    else
    {
      // If dimension is fixed
      if(E_.norm() == 0)
      {
        Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZero;
        ABZero << dt_ * A_, dt_ * B_, Eigen::Matrix<double, InputDim, StateDim + InputDim>::Zero();
        Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZeroExp = ABZero.exp();
        Ad_ = ABZeroExp.template block<StateDim, StateDim>(0, 0);
        Bd_ = ABZeroExp.template block<StateDim, InputDim>(0, StateDim);
        Ed_.setZero();
      }
      else
      {
        // There is no proof that this is correct.
        Eigen::Matrix<double, StateDim + InputDim + 1, StateDim + InputDim + 1> ABEZero;
        ABEZero << dt_ * A_, dt_ * B_, dt_ * E_, Eigen::Matrix<double, InputDim + 1, StateDim + InputDim + 1>::Zero();
        Eigen::Matrix<double, StateDim + InputDim + 1, StateDim + InputDim + 1> ABEZeroExp = ABEZero.exp();
        Ad_ = ABEZeroExp.template block<StateDim, StateDim>(0, 0);
        Bd_ = ABEZeroExp.template block<StateDim, InputDim>(0, StateDim);
        Ed_ = ABEZeroExp.template block<StateDim, 1>(0, StateDim + InputDim);
      }
    }
  }

public:
  //! State dimension
  const int state_dim_ = 0;

  //! Input dimension
  const int input_dim_ = 0;

  //! Output dimension
  const int output_dim_ = 0;

  //! Matrix A of continuous state equation
  Eigen::Matrix<double, StateDim, StateDim> A_;

  //! Matrix B of continuous state equation
  Eigen::Matrix<double, StateDim, InputDim> B_;

  //! Matrix C of observation equation
  Eigen::Matrix<double, OutputDim, StateDim> C_;

  //! Matrix D of observation equation
  Eigen::Matrix<double, OutputDim, InputDim> D_;

  //! Offset vector E of continuous state equation
  StateDimVector E_;

  //! Offset vector F of observation equation
  OutputDimVector F_;

  //! Discretization timestep [sec] (zero if discrete coefficients are not initialized)
  double dt_ = 0;

  //! Matrix A of discrete state equation
  Eigen::Matrix<double, StateDim, StateDim> Ad_;

  //! Matrix B of discrete state equation
  Eigen::Matrix<double, StateDim, InputDim> Bd_;

  //! Offset vector E of discrete state equation
  StateDimVector Ed_;
};
} // namespace CCC
