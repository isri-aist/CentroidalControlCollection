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
    if constexpr (StateDim == Eigen::Dynamic) {
        return A_.rows();
      } else {
      return StateDim;
    }
  }

  /*! \brief Gets the input dimension. */
  int inputDim() const
  {
    if constexpr (InputDim == Eigen::Dynamic) {
        return B_.cols();
      } else {
      return InputDim;
    }
  }

  /*! \brief Gets the output dimension. */
  int outputDim() const
  {
    if constexpr (OutputDim == Eigen::Dynamic) {
        return C_.rows();
      } else {
      return OutputDim;
    }
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
    if constexpr (StateDim == Eigen::Dynamic || InputDim == Eigen::Dynamic) {
      if (E_.norm() == 0) {
        Eigen::MatrixXd ABZero(stateDim() + inputDim(), stateDim() + inputDim());
        ABZero << dt_ * A_, dt_ * B_, Eigen::MatrixXd::Zero(inputDim(), stateDim() + inputDim());
        Eigen::MatrixXd ABZeroExp = ABZero.exp();
        Ad_ = ABZeroExp.template block(0, 0, stateDim(), stateDim());
        Bd_ = ABZeroExp.template block(0, stateDim(), stateDim(), inputDim());
      } else {
        // There is no proof that this is correct.
        Eigen::Matrix<double, StateDim + InputDim + 1, StateDim + InputDim + 1> ABEZero;
        ABEZero << dt_ * A_, dt_ * B_, dt_ * E_, Eigen::Matrix<double, InputDim + 1, StateDim + InputDim + 1>::Zero();
        Eigen::Matrix<double, StateDim + InputDim + 1, StateDim + InputDim + 1> ABEZeroExp = ABEZero.exp();
        Ad_ = ABEZeroExp.template block(0, 0, stateDim(), stateDim());
        Bd_ = ABEZeroExp.template block(0, stateDim(), stateDim(), inputDim());
        Ed_ = ABEZeroExp.template block(0, stateDim() + inputDim(), stateDim(), 1);
      }
      } else {
      if (E_.norm() == 0) {
        Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZero;
        ABZero << dt_ * A_, dt_ * B_, Eigen::Matrix<double, InputDim, StateDim + InputDim>::Zero();
        Eigen::Matrix<double, StateDim + InputDim, StateDim + InputDim> ABZeroExp = ABZero.exp();
        Ad_ = ABZeroExp.template block<StateDim, StateDim>(0, 0);
        Bd_ = ABZeroExp.template block<StateDim, InputDim>(0, StateDim);
      } else {
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
  //! Matrix A of continuous state equation
  Eigen::Matrix<double, StateDim, StateDim> A_;

  //! Matrix B of continuous state equation
  Eigen::Matrix<double, StateDim, InputDim> B_;

  //! Matrix C of observation equation
  Eigen::Matrix<double, OutputDim, StateDim> C_;

  //! Matrix D of observation equation
  Eigen::Matrix<double, OutputDim, InputDim> D_;

  //! Offset vector of continuous state equation
  StateDimVector E_;

  //! Timestep for discretization (-1 if discretization coefficients are not initialized)
  double dt_ = -1;

  //! Matrix A of discrete state equation
  Eigen::Matrix<double, StateDim, StateDim> Ad_;

  //! Matrix B of discrete state equation
  Eigen::Matrix<double, StateDim, InputDim> Bd_;

  //! Offset vector of discrete state equation
  StateDimVector Ed_;
};
}
