/* Author: Masaki Murooka */

#pragma once

#include <memory>

#include <CCC/StateSpaceModel.h>
#include <CCC/console.h>

namespace CCC
{
/** \brief Preview control.
    \tparam StateDim state dimension
    \tparam InputDim input dimension
    \tparam OutputDim output dimension

    See the following for a detailed formulation.
      - Shuuji Kajita, et al. Biped walking pattern generation by using preview control of zero-moment point. ICRA,
   2003.
 */
template<int StateDim, int InputDim, int OutputDim>
class PreviewControl
{
public:
  /** \brief Type of state vector. */
  using StateDimVector = Eigen::Matrix<double, StateDim, 1>;

  /** \brief Type of input vector. */
  using InputDimVector = Eigen::Matrix<double, InputDim, 1>;

  /** \brief Type of output vector. */
  using OutputDimVector = Eigen::Matrix<double, OutputDim, 1>;

  /** \brief Type of state-space model. */
  using _StateSpaceModel = StateSpaceModel<StateDim, InputDim, OutputDim>;

public:
  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! Output weight
    OutputDimVector output;

    //! Input weight
    InputDimVector input;

    /** \brief Constructor.
        \param _output output weight
        \param _input input weight
     */
    WeightParam(const OutputDimVector & _output = OutputDimVector::Zero(),
                const InputDimVector & _input = InputDimVector::Zero())
    : output(_output), input(_input)
    {
    }
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param model state-space model
      \param horizon_duration horizon duration [sec]
      \param horizon_dt discretization timestep in horizon [sec]
      \param weight_param objective weight parameter
   */
  PreviewControl(const std::shared_ptr<_StateSpaceModel> & model,
                 double horizon_duration,
                 double horizon_dt,
                 const WeightParam & weight_param)
  : model_(model), horizon_dt_(horizon_dt), horizon_steps_(static_cast<int>(std::ceil(horizon_duration / horizon_dt)))
  {
    if(horizon_duration <= 0 || horizon_dt <= 0)
    {
      throw std::runtime_error("[PreviewControl] Input arguments are invalid. horizon_duration: "
                               + std::to_string(horizon_duration) + ", horizon_dt: " + std::to_string(horizon_dt));
    }

    calcGain(weight_param);
  }

  /** \brief Calculate optimal input.
      \param x state
      \param ref_output_seq sequence of reference output
   */
  InputDimVector calcOptimalInput(const StateDimVector & x, const Eigen::VectorXd & ref_output_seq) const
  {
    return -K_ * x + F_ * ref_output_seq;
  }

protected:
  /** \brief Calculate the gain of the preview control. */
  void calcGain(const WeightParam & weight_param)
  {
    // 0. Calculate discrete system matrices
    if(model_->dt_ != horizon_dt_)
    {
      model_->calcDiscMatrix(horizon_dt_);
    }

    const Eigen::Matrix<double, StateDim, StateDim> & A = model_->Ad_;
    const Eigen::Matrix<double, StateDim, InputDim> & B = model_->Bd_;
    const Eigen::Matrix<double, OutputDim, StateDim> & C = model_->C_;

    Eigen::Matrix<double, OutputDim, OutputDim> Q = weight_param.output.asDiagonal();
    Eigen::Matrix<double, InputDim, InputDim> R = weight_param.input.asDiagonal();
    Eigen::Matrix<double, InputDim, InputDim> Rinv = weight_param.input.cwiseInverse().asDiagonal();

    // int stateDim = model_->stateDim();
    int inputDim = model_->inputDim();
    int outputDim = model_->outputDim();

    // 1. Calculate P by solving the Discrete-time Algebraic Riccati Equation
    // (DARE) ref. https://scicomp.stackexchange.com/a/35394
    int iterMax = 10000;
    double relNormThre = 1e-8;
    Eigen::Matrix<double, StateDim, StateDim> A0 = A;
    Eigen::Matrix<double, StateDim, StateDim> G0 = B * Rinv * B.transpose();
    Eigen::Matrix<double, StateDim, StateDim> H0 = C.transpose() * Q * C;
    Eigen::Matrix<double, StateDim, StateDim> A1, G1, H1;
    for(int i = 0; i < iterMax; i++)
    {
      Eigen::Matrix<double, StateDim, StateDim> I_GH_inv =
          (Eigen::Matrix<double, StateDim, StateDim>::Identity() + G0 * H0).inverse();
      A1 = A0 * I_GH_inv * A0;
      G1 = G0 + A0 * I_GH_inv * G0 * A0.transpose();
      H1 = H0 + A0.transpose() * H0 * I_GH_inv * A0;

      double relNorm = (H1 - H0).norm() / H1.norm();
      if(relNorm < relNormThre)
      {
        break;
      }
      else if(i == iterMax - 1)
      {
        CCC_WARN_STREAM("[PreviewControl] Solution of Riccati equation did not converged: " << relNorm << " > "
                                                                                            << relNormThre);
      }

      A0 = A1;
      G0 = G1;
      H0 = H1;
    }
    P_ = H1;

    riccati_error_ = (P_
                      - (A.transpose() * P_ * A + C.transpose() * Q * C
                         - A.transpose() * P_ * B * (R + B.transpose() * P_ * B).inverse() * B.transpose() * P_ * A))
                         .norm();

    // 2. Calculate the feedback gain (K and f)
    Eigen::Matrix<double, InputDim, InputDim> R_BtPB_inv = (R + B.transpose() * P_ * B).inverse();
    // 2.1 Calculate K
    K_ = R_BtPB_inv * B.transpose() * P_ * A;
    // 2.2 Calculate f
    F_.resize(inputDim, horizon_steps_ * outputDim);
    Eigen::Matrix<double, StateDim, StateDim> A_BK = A - B * K_;
    Eigen::Matrix<double, StateDim, StateDim> fSub = Eigen::Matrix<double, StateDim, StateDim>::Identity();
    for(int i = 0; i < horizon_steps_; i++)
    {
      if(i < horizon_steps_ - 1)
      {
        F_.middleCols(i * outputDim, outputDim) = R_BtPB_inv * B.transpose() * fSub * C.transpose() * Q;
      }
      else
      {
        // See https://github.com/euslisp/jskeus/pull/551
        F_.middleCols(i * outputDim, outputDim) = R_BtPB_inv * B.transpose() * fSub * P_ * C.transpose();
      }
      fSub = fSub * A_BK.transpose();
    }
  }

public:
  //! State-space model
  std::shared_ptr<_StateSpaceModel> model_;

  //! Discretization timestep in horizon [sec]
  double horizon_dt_ = 0;

  //! Number of steps in horizon
  int horizon_steps_ = -1;

  //! Error of algebraic Riccati equation
  double riccati_error_ = 0;

protected:
  //! Solution of the algebraic Riccati equation
  Eigen::Matrix<double, StateDim, StateDim> P_;

  //! Feedback gain
  Eigen::Matrix<double, InputDim, StateDim> K_;

  //! Preview gain
  Eigen::Matrix<double, InputDim, Eigen::Dynamic> F_;
};
} // namespace CCC
