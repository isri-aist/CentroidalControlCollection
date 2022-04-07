/* Author: Masaki Murooka */

#pragma once

#include <memory>

#include <ros/console.h>

#include <CCC/StateSpaceModel.h>

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
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param model state-space model
   */
  PreviewControl(const std::shared_ptr<_StateSpaceModel> & model) : model_(model) {}

  /** \brief Calculate the gain of the preview control.
      \param horizon horizon of the preview control [sec]
      \param dt sampling time of the preview control [sec]
      \param output_weight weight of output
      \param input_weight weight of input
   */
  void calcGain(double horizon, double dt, OutputDimVector output_weight, InputDimVector input_weight)
  {
    if(horizon <= 0 || dt <= 0)
    {
      throw std::runtime_error("Input arguments of PreviewControl::calcGain are invalid. horizon: "
                               + std::to_string(horizon) + ", dt: " + std::to_string(dt));
    }

    // 0. Calculate discrete system matrices
    if(model_->dt_ != dt)
    {
      model_->calcDiscMatrix(dt);
    }

    const Eigen::Matrix<double, StateDim, StateDim> & A = model_->Ad_;
    const Eigen::Matrix<double, StateDim, InputDim> & B = model_->Bd_;
    const Eigen::Matrix<double, OutputDim, StateDim> & C = model_->C_;

    Eigen::Matrix<double, OutputDim, OutputDim> Q = output_weight.asDiagonal();
    Eigen::Matrix<double, InputDim, InputDim> R = input_weight.asDiagonal();
    Eigen::Matrix<double, InputDim, InputDim> Rinv = input_weight.cwiseInverse().asDiagonal();

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
        ROS_WARN_STREAM("[PreviewControl] Solution of Riccati equation did not converged: " << relNorm << " > "
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
    preview_size_ = static_cast<int>(std::ceil(horizon / dt));
    Eigen::Matrix<double, InputDim, InputDim> R_BtPB_inv = (R + B.transpose() * P_ * B).inverse();
    // 2.1 Calculate K
    K_ = R_BtPB_inv * B.transpose() * P_ * A;
    // 2.2 Calculate f
    F_.resize(inputDim, preview_size_ * outputDim);
    Eigen::Matrix<double, StateDim, StateDim> A_BK = A - B * K_;
    Eigen::Matrix<double, StateDim, StateDim> fSub = Eigen::Matrix<double, StateDim, StateDim>::Identity();
    for(int i = 0; i < preview_size_; i++)
    {
      if(i < preview_size_ - 1)
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

  /** \brief Calculate optimal input.
      \param x state
      \param ref_output_traj trajectory of reference output
   */
  InputDimVector calcOptimalInput(const StateDimVector & x, const Eigen::VectorXd & ref_output_traj) const
  {
    return -K_ * x + F_ * ref_output_traj;
  }

public:
  //! State-space model
  std::shared_ptr<_StateSpaceModel> model_;

  //! Solution of the algebraic Riccati equation
  Eigen::Matrix<double, StateDim, StateDim> P_;

  //! Feedback gain
  Eigen::Matrix<double, InputDim, StateDim> K_;

  //! Preview gain
  Eigen::Matrix<double, InputDim, Eigen::Dynamic> F_;

  //! Error of algebraic Riccati equation
  double riccati_error_ = 0;

  //! Size of preview window
  int preview_size_ = -1;
};
} // namespace CCC
