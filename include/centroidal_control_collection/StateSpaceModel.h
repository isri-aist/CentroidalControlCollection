#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace CCC
{
/*! \brief State-space model.
 *
 *  \tparam StateDim state dimension
 *  \tparam InputDim input dimension
 *  \tparam OutputDim output dimension
 */
template<int StateDim, int InputDim, int OutputDim>
class StateSpaceModel
{
public:
  /*! \brief Type of state vector. */
  using StateDimVector = Eigen::Matrix<double, StateDim, 1>;

  /*! \brief Type of input vector. */
  using InputDimVector = Eigen::Matrix<double, InputDim, 1>;

  /*! \brief Type of output vector. */
  using OutputDimVector = Eigen::Matrix<double, OutputDim, 1>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*! \brief Constructor. */
  StateSpaceModel() {}

  /*! \brief Destructor.
   *
   *  \note Need to make class polymorphic.
   *        See https://stackoverflow.com/a/15114118
   */
  virtual ~StateSpaceModel() = default;

  /*! \brief Gets the state dimension. */
  constexpr int stateDim() const
  {
    return StateDim;
  }

  /*! \brief Gets the input dimension. */
  constexpr int inputDim() const
  {
    return InputDim;
  }

  /*! \brief Gets the output dimension. */
  constexpr int outputDim() const
  {
    return OutputDim;
  }

  /*! \brief Calculate the result of the continuous state equation.
   *
   *  \param x state
   *  \param u input
   *  \returns time derivative of x (\dot{x})
   */
  StateDimVector stateEq(const StateDimVector & x, const InputDimVector & u) const
  {
    return A_ * x + B_ * u; // dx
  }

  /*! \brief Calculate the result of the discrete state equation.
   *
   *  \param x current state (x[k])
   *  \param u current input (u[k])
   *  \returns next state (x[k+1])
   */
  StateDimVector stateEqDisc(const StateDimVector & x, const InputDimVector & u) const
  {
    return Ad_ * x + Bd_ * u; // next_x
  }

  /*! \brief Calculate the result of the observation equation.
   *
   *  \param x state
   */
  OutputDimVector observEq(const StateDimVector & x) const
  {
    return C_ * x; // y
  }

  /*! \brief Calculate the result of the observation equation.
   *
   *  \param x state
   *  \param u input
   */
  OutputDimVector observEq(const StateDimVector & x, const InputDimVector & u) const
  {
    return C_ * x + D_ * u; // y
  }

  /*! \brief Calculate the discrete system matrices.
   *
   *  \param dt timestep
   */
  void calcDiscMatrix(double dt)
  {
    dt_ = dt;

    // Zero-order hold discretization
    // See https://en.wikipedia.org/wiki/Discretization
    Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZero;
    ABZero << dt_ * A_, dt_ * B_, Eigen::Matrix<double, InputDim, StateDim + InputDim>::Zero();
    Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZeroExp = ABZero.exp();
    Ad_ = ABZeroExp.block(0, 0, A_.rows(), A_.cols());
    Bd_ = ABZeroExp.block(0, A_.cols(), B_.rows(), B_.cols());
  }

public:
  //! Matrix A of continuous state equation
  Eigen::Matrix<double, StateDim, StateDim> A_;

  //! Matrix B of continuous state equation
  Eigen::Matrix<double, StateDim, InputDim> B_;

  //! Matrix C of observation equation
  Eigen::Matrix<double, OutputDim, StateDim> C_;

  //! Matrix D of observation equation
  Eigen::Matrix<double, OutputDim, InputDim> D_;

  //! Timestep for discretization
  double dt_ = -1;

  //! Matrix A of discrete state equation
  Eigen::Matrix<double, StateDim, StateDim> Ad_;

  //! Matrix B of discrete state equation
  Eigen::Matrix<double, StateDim, InputDim> Bd_;
};
} // namespace Motion6DoF