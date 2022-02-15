/* Author: Masaki Murooka */

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
  int stateDim() const
  {
    return StateDim;
  }

  /*! \brief Gets the input dimension. */
  int inputDim() const
  {
    return InputDim;
  }

  /*! \brief Gets the output dimension. */
  int outputDim() const
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
    return A_ * x + B_ * u + E_; // dx
  }

  /*! \brief Calculate the result of the discrete state equation.
   *
   *  \param x current state (x[k])
   *  \param u current input (u[k])
   *  \returns next state (x[k+1])
   */
  StateDimVector stateEqDisc(const StateDimVector & x, const InputDimVector & u) const
  {
    return Ad_ * x + Bd_ * u + Ed_; // next_x
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
    if (E_.norm() == 0) {
      Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZero;
      ABZero << dt_ * A_, dt_ * B_, Eigen::Matrix<double, InputDim, StateDim + InputDim>::Zero();
      Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZeroExp = ABZero.exp();
      Ad_ = ABZeroExp.block(0, 0, A_.rows(), A_.cols());
      Bd_ = ABZeroExp.block(0, A_.cols(), B_.rows(), B_.cols());
    } else {
      // There is no proof that this is correct.
      Eigen::Matrix<double, StateDim + InputDim + 1, StateDim + InputDim + 1> ABEZero;
      ABEZero << dt_ * A_, dt_ * B_, dt_ * E_, Eigen::Matrix<double, InputDim + 1, StateDim + InputDim + 1>::Zero();
      Eigen::Matrix<double, StateDim + InputDim + 1, StateDim + InputDim + 1> ABEZeroExp = ABEZero.exp();
      Ad_ = ABEZeroExp.block(0, 0, A_.rows(), A_.cols());
      Bd_ = ABEZeroExp.block(0, A_.cols(), B_.rows(), B_.cols());
      Ed_ = ABEZeroExp.block(0, A_.cols() + B_.cols(), E_.rows(), E_.cols());
    }
  }

public:
  //! Matrix A of continuous state equation
  Eigen::Matrix<double, StateDim, StateDim> A_ = Eigen::Matrix<double, StateDim, StateDim>::Zero();

  //! Matrix B of continuous state equation
  Eigen::Matrix<double, StateDim, InputDim> B_ = Eigen::Matrix<double, StateDim, InputDim>::Zero();

  //! Offset vector of continuous state equation
  StateDimVector E_ = StateDimVector::Zero();

  //! Matrix C of observation equation
  Eigen::Matrix<double, OutputDim, StateDim> C_ = Eigen::Matrix<double, OutputDim, StateDim>::Zero;

  //! Matrix D of observation equation
  Eigen::Matrix<double, OutputDim, InputDim> D_ = Eigen::Matrix<double, OutputDim, InputDim>::Zero();

  //! Timestep for discretization (-1 if discretization coefficients are not initialized)
  double dt_ = -1;

  //! Matrix A of discrete state equation
  Eigen::Matrix<double, StateDim, StateDim> Ad_ = Eigen::Matrix<double, StateDim, StateDim>::Zero();

  //! Matrix B of discrete state equation
  Eigen::Matrix<double, StateDim, InputDim> Bd_ = Eigen::Matrix<double, StateDim, InputDim>::Zero();

  //! Offset vector of discrete state equation
  StateDimVector Ed_ = StateDimVector::Zero();
};

template <>
inline int StateSpaceModel<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic>::stateDim() const
{
  return A_.rows();
}

template <>
inline int StateSpaceModel<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic>::inputDim() const
{
  return B_.cols();
}

template <>
inline int StateSpaceModel<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic>::outputDim() const
{
  return C_.rows();
}

using StateSpaceModelDynamicShape = StateSpaceModel<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic>;
}
