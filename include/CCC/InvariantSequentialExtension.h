/* Author: Masaki Murooka */

#pragma once

#include <iostream>
#include <memory>

#include <CCC/StateSpaceModel.h>

namespace CCC
{
/** \brief Sequential extension for time-invariant system.
    \tparam StateDim state dimension
    \tparam InputDim input dimension
    \tparam OutputDim output dimension

    Given the following time-invariant linear discrete state equation.
    \f{align*}{ \boldsymbol{x}_{k+1} = \boldsymbol{A} \boldsymbol{x}_{k} +
    \boldsymbol{B} \boldsymbol{u}_{k} + \boldsymbol{e} \f}

    The following equation is called "sequential extension" here. In this class, the coefficients
   \f$\boldsymbol{\hat{A}}\f$, \f$\boldsymbol{\hat{B}}\f$, and \f$\boldsymbol{\hat{e}}\f$ are calculated.
    \f{align*}{
    \boldsymbol{\hat{x}}_{k+1} &= \boldsymbol{\hat{A}} \boldsymbol{x}_{k} + \boldsymbol{\hat{B}}
   \boldsymbol{\hat{u}}_{k} + \boldsymbol{\hat{e}} \\ \Leftrightarrow \begin{bmatrix} \boldsymbol{x}_{k+1} \\
   \boldsymbol{x}_{k+2} \\ \boldsymbol{x}_{k+3} \\ \vdots \\ \boldsymbol{x}_{k+N} \end{bmatrix}
   &= \begin{bmatrix}
      \boldsymbol{A} \\
      \boldsymbol{A}^2 \\
      \boldsymbol{A}^3 \\
      \vdots \\
      \boldsymbol{A}^N
    \end{bmatrix}
    \boldsymbol{x}_{k} +
    \begin{bmatrix}
      \boldsymbol{B} & \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} & \cdots & \boldsymbol{O} \\
      \boldsymbol{A} \boldsymbol{B} & \boldsymbol{B} & \boldsymbol{O} & \boldsymbol{O} & \cdots & \boldsymbol{O} \\
      \boldsymbol{A}^2 \boldsymbol{B} & \boldsymbol{A} \boldsymbol{B} & \boldsymbol{B} & \boldsymbol{O} & \cdots &
   \boldsymbol{O} \\
      \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
      \boldsymbol{A}^{N-1} \boldsymbol{B} & \cdots & \cdots & \cdots & \boldsymbol{A} \boldsymbol{B} & \boldsymbol{B}
    \end{bmatrix} \begin{bmatrix}
   \boldsymbol{u}_{k}
   \\ \boldsymbol{u}_{k+1} \\ \boldsymbol{u}_{k+2} \\ \vdots \\ \boldsymbol{u}_{k+N-1}
    \end{bmatrix} +
    \begin{bmatrix}
      \boldsymbol{e} \\
      (\boldsymbol{A} + \boldsymbol{I}) \boldsymbol{e} \\
      (\boldsymbol{A}^2 + \boldsymbol{A} + \boldsymbol{I}) \boldsymbol{e} \\
      \vdots \\
      (\boldsymbol{A}^{N-1} + \cdots + \boldsymbol{A} + \boldsymbol{I}) \boldsymbol{e}
   \end{bmatrix}
   \f}

    Such a sequential extension is often used to formulate linear MPC as quadratic programming. For example, the
   following papers uses it.
      - PB Wieber. Trajectory Free Linear Model Predictive Control for Stable Walking in the Presence of Strong
   Perturbations. Humanoids, 2006.
*/
template<int StateDim, int InputDim, int OutputDim>
class InvariantSequentialExtension
{
public:
  /** \brief Type of state-space model with fixed dimensions of state, input, and output. */
  using _StateSpaceModel = StateSpaceModel<StateDim, InputDim, OutputDim>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param model state-space model
      \param seq_len sequence length
      \param extend_for_output whether to extend for output instead of state
   */
  InvariantSequentialExtension(const std::shared_ptr<_StateSpaceModel> & model,
                               int seq_len,
                               bool extend_for_output = false)
  : model_(model), seq_len_(seq_len)
  {
    setup(extend_for_output);
  }

  /** \brief Get total state dimension. */
  int totalStateDim() const
  {
    return seq_len_ * StateDim;
  }

  /** \brief Get total input dimension. */
  int totalInputDim() const
  {
    return seq_len_ * InputDim;
  }

  /** \brief Get total output dimension. */
  int totalOutputDim() const
  {
    return seq_len_ * OutputDim;
  }

protected:
  /** \brief Setup coefficients. */
  void setup(bool extend_for_output)
  {
    // Check D matrix
    if(extend_for_output && model_->D_.norm() > 0)
    {
      std::cerr
          << "[InvariantSequentialExtension] Matrix D in state-space model must be zero when extending for output.\n"
          << "  Matrix D:\n"
          << model_->D_;
    }

    // Resize matrix
    A_seq_.setZero(seq_len_ * StateDim, StateDim);
    B_seq_.setZero(seq_len_ * StateDim, seq_len_ * InputDim);
    E_seq_.setZero(seq_len_ * StateDim);

    for(int i = 0; i < seq_len_; i++)
    {
      // Setq A_seq_
      if(i == 0)
      {
        A_seq_.template middleRows<StateDim>(i * StateDim) = model_->Ad_;
      }
      else
      {
        (A_seq_.template middleRows<StateDim>(i * StateDim)).noalias() =
            model_->Ad_ * A_seq_.template middleRows<StateDim>((i - 1) * StateDim);
      }

      // Setq B_seq_
      for(int j = 0; j < seq_len_ - i; j++)
      {
        if(j == 0)
        {
          if(i == 0)
          {
            B_seq_.template block<StateDim, InputDim>(i * StateDim, 0) = model_->Bd_;
          }
          else
          {
            B_seq_.template block<StateDim, InputDim>(i * StateDim, 0).noalias() =
                model_->Ad_ * B_seq_.template block<StateDim, InputDim>((i - 1) * StateDim, 0);
          }
        }
        else
        {
          B_seq_.template block<StateDim, InputDim>((i + j) * StateDim, j * InputDim) =
              B_seq_.template block<StateDim, InputDim>(i * StateDim, 0);
        }
      }

      // Setq E_seq_
      if(i == 0)
      {
        E_seq_.template segment<StateDim>(i * StateDim) = model_->Ed_;
      }
      else
      {
        (E_seq_.template segment<StateDim>(i * StateDim)).noalias() =
            model_->Ad_ * E_seq_.template segment<StateDim>((i - 1) * StateDim) + model_->Ed_;
      }
    }

    // Transform for output
    if(extend_for_output)
    {
      // Set C_seq
      Eigen::MatrixXd C_seq = Eigen::MatrixXd::Zero(seq_len_ * OutputDim, seq_len_ * StateDim);
      for(int i = 0; i < seq_len_; i++)
      {
        C_seq.template block<OutputDim, StateDim>(i * OutputDim, i * StateDim) = model_->C_;
      }

      // Apply C_seq
      A_seq_ = C_seq * A_seq_;
      B_seq_ = C_seq * B_seq_;
      E_seq_ = C_seq * E_seq_;
    }
  }

public:
  //! State-space model
  std::shared_ptr<_StateSpaceModel> model_;

  //! Sequence length
  int seq_len_ = 0;

  //! Sequential extension matrix \f$\boldsymbol{\hat{A}}\f$
  Eigen::Matrix<double, Eigen::Dynamic, StateDim> A_seq_;

  //! Sequential extension matrix \f$\boldsymbol{\hat{B}}\f$
  Eigen::MatrixXd B_seq_;

  //! Sequential extension vector (i.e., offset vector) \f$\boldsymbol{\hat{e}}\f$
  Eigen::VectorXd E_seq_;
};
} // namespace CCC
