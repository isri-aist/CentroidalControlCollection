/* Author: Masaki Murooka */

#pragma once

#include <iostream>
#include <memory>

#include <CCC/StateSpaceModel.h>

namespace CCC
{
/** \brief Sequential extension for time-variant system.
    \tparam StateDim state dimension
    \tparam ListType type of state-space model list
    \note State dimension must be the same for all models in the sequence.

    Given the following time-variant linear discrete state equation (even systems with time-variant dimensions of
   control input are acceptable). \f{align*}{ \boldsymbol{x}_{k+1} = \boldsymbol{A}_{k} \boldsymbol{x}_{k} +
   \boldsymbol{B}_{k} \boldsymbol{u}_{k} + \boldsymbol{e}_k \f}

    The following equation is called "sequential extension" here. In this class, the coefficients
   \f$\boldsymbol{\hat{A}}_k\f$, \f$\boldsymbol{\hat{B}}_k\f$, and \f$\boldsymbol{\hat{e}}_k\f$ are calculated.
    \f{align*}{
    \boldsymbol{\hat{x}}_{k+1} &= \boldsymbol{\hat{A}}_{k} \boldsymbol{x}_{k} + \boldsymbol{\hat{B}}_{k}
   \boldsymbol{\hat{u}}_{k} + \boldsymbol{\hat{e}}_k \\ \Leftrightarrow \begin{bmatrix} \boldsymbol{x}_{k+1} \\
   \boldsymbol{x}_{k+2} \\ \boldsymbol{x}_{k+3} \\ \boldsymbol{x}_{k+4} \\ \vdots \\ \boldsymbol{x}_{k+N} \end{bmatrix}
   &= \begin{bmatrix}
      \boldsymbol{A}_{k} \\
      \boldsymbol{A}_{k+1} \boldsymbol{A}_{k} \\
      \boldsymbol{A}_{k+2} \boldsymbol{A}_{k+1} \boldsymbol{A}_{k} \\
      \boldsymbol{A}_{k+3} \boldsymbol{A}_{k+2} \boldsymbol{A}_{k+1} \boldsymbol{A}_{k} \\
      \vdots \\
      \boldsymbol{A}_{k+N-1} \cdots \boldsymbol{A}_{k+1} \boldsymbol{A}_{k}
    \end{bmatrix}
    \boldsymbol{x}_{k} \\
    & +
    \begin{bmatrix}
      \boldsymbol{B}_{k} & \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} & \cdots & \boldsymbol{O}
   \\
      \boldsymbol{A}_{k+1} \boldsymbol{B}_{k} & \boldsymbol{B}_{k+1} & \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O}
   & \cdots & \boldsymbol{O} \\
      \boldsymbol{A}_{k+2} \boldsymbol{A}_{k+1} \boldsymbol{B}_{k} & \boldsymbol{A}_{k+2} \boldsymbol{B}_{k+1} &
   \boldsymbol{B}_{k+2} & \boldsymbol{O} & \boldsymbol{O} & \cdots & \boldsymbol{O} \\
      \boldsymbol{A}_{k+3} \boldsymbol{A}_{k+2} \boldsymbol{A}_{k+1} \boldsymbol{B}_{k} & \boldsymbol{A}_{k+3}
   \boldsymbol{A}_{k+2} \boldsymbol{B}_{k+1} & \boldsymbol{A}_{k+3} \boldsymbol{B}_{k+2} & \boldsymbol{B}_{k+3} &
   \boldsymbol{O} & \cdots & \boldsymbol{O} \\
      \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
      \boldsymbol{A}_{k+N-1} \cdots \boldsymbol{A}_{k+1} \boldsymbol{B}_{k} & \cdots & \cdots & \cdots & \cdots &
   \boldsymbol{A}_{k+N-1} \boldsymbol{B}_{k+N-2} & \boldsymbol{B}_{k+N-1} \end{bmatrix} \begin{bmatrix}
   \boldsymbol{u}_{k}
   \\ \boldsymbol{u}_{k+1} \\ \boldsymbol{u}_{k+2} \\ \boldsymbol{u}_{k+3} \\ \vdots \\ \boldsymbol{u}_{k+N-1}
    \end{bmatrix} \\
    & +
    \begin{bmatrix}
      \boldsymbol{e}_{k} \\
      \boldsymbol{A}_{k+1} \boldsymbol{e}_{k} + \boldsymbol{e}_{k+1} \\
      \boldsymbol{A}_{k+2} \boldsymbol{A}_{k+1} \boldsymbol{e}_{k} + \boldsymbol{A}_{k+2} \boldsymbol{e}_{k+1} +
   \boldsymbol{e}_{k+2}\\
      \vdots \\
      \vdots \\
      \boldsymbol{A}_{k+N-1} \cdots \boldsymbol{A}_{k+1} \boldsymbol{e}_{k} +
   \boldsymbol{A}_{k+N-1} \cdots\boldsymbol{A}_{k+2} \boldsymbol{e}_{k+1} + \cdots +
   \boldsymbol{e}_{k+N-1} \end{bmatrix} \f}

    Such a sequential extension is often used to formulate linear MPC as quadratic programming. For example, the
   following papers uses it.
      - H Audren, et al. Model preview control in multi-contact motion-application to a humanoid robot. IROS, 2014.
*/
template<int StateDim, template<class> class ListType = std::vector>
class VariantSequentialExtension
{
public:
  /** \brief Type of state-space model with fixed state dimension. */
  using _StateSpaceModel = StateSpaceModel<StateDim, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** \brief Constructor.
      \param model_list list of state-space model
      \param extend_for_output whether to extend for output instead of state
   */
  VariantSequentialExtension(const ListType<std::shared_ptr<_StateSpaceModel>> & model_list,
                             bool extend_for_output = false)
  : model_list_(model_list)
  {
    setup(extend_for_output);
  }

  /** \brief Get total state dimension. */
  int totalStateDim() const
  {
    return total_state_dim_;
  }

  /** \brief Get total input dimension. */
  int totalInputDim() const
  {
    return total_input_dim_;
  }

  /** \brief Get total output dimension. */
  int totalOutputDim() const
  {
    return total_output_dim_;
  }

protected:
  /** \brief Setup coefficients. */
  void setup(bool extend_for_output)
  {
    // Check D matrix
    if(extend_for_output)
    {
      for(const auto & model : model_list_)
      {
        if(model->D_.norm() > 0)
        {
          std::cerr
              << "[VariantSequentialExtension] Matrix D in state-space model must be zero when extending for output.\n"
              << "  Matrix D:\n"
              << model->D_;
          break;
        }
      }
    }

    // Set total dimension
    size_t seq_len = model_list_.size();
    total_state_dim_ = seq_len * StateDim;
    total_input_dim_ = 0;
    total_output_dim_ = 0;
    for(const auto & model : model_list_)
    {
      total_input_dim_ += model->inputDim();
      total_output_dim_ += model->outputDim();
    }

    // Resize matrix
    A_seq_.setZero(total_state_dim_, StateDim);
    B_seq_.setZero(total_state_dim_, total_input_dim_);
    E_seq_.setZero(total_state_dim_);

    int accum_input_dim = 0;
    for(size_t i = 0; i < seq_len; i++)
    {
      int current_input_dim = model_list_[i]->inputDim();

      // Setq A_seq_
      if(i == 0)
      {
        A_seq_.template middleRows<StateDim>(i * StateDim) = model_list_[0]->Ad_;
      }
      else
      {
        (A_seq_.template middleRows<StateDim>(i * StateDim)).noalias() =
            model_list_[i]->Ad_ * A_seq_.template middleRows<StateDim>((i - 1) * StateDim);
      }

      // Setq B_seq_
      for(size_t j = i; j < seq_len; j++)
      {
        if(j == i)
        {
          B_seq_.block(j * StateDim, accum_input_dim, StateDim, current_input_dim) = model_list_[i]->Bd_;
        }
        else
        {
          (B_seq_.block(j * StateDim, accum_input_dim, StateDim, current_input_dim)).noalias() =
              model_list_[j]->Ad_ * B_seq_.block((j - 1) * StateDim, accum_input_dim, StateDim, current_input_dim);
        }
      }

      // Setq E_seq_
      if(i == 0)
      {
        E_seq_.segment<StateDim>(i * StateDim) = model_list_[0]->Ed_;
      }
      else
      {
        (E_seq_.segment<StateDim>(i * StateDim)).noalias() =
            model_list_[i]->Ad_ * E_seq_.segment<StateDim>((i - 1) * StateDim) + model_list_[i]->Ed_;
      }

      accum_input_dim += current_input_dim;
    }

    // Transform for output
    if(extend_for_output)
    {
      // Set C_seq
      Eigen::MatrixXd C_seq = Eigen::MatrixXd::Zero(total_output_dim_, total_state_dim_);
      int accum_output_dim = 0;
      for(size_t i = 0; i < seq_len; i++)
      {
        int current_output_dim = model_list_[i]->outputDim();

        C_seq.block(accum_output_dim, i * StateDim, current_output_dim, StateDim) = model_list_[i]->C_;

        accum_output_dim += current_output_dim;
      }

      // Apply C_seq
      A_seq_ = C_seq * A_seq_;
      B_seq_ = C_seq * B_seq_;
      E_seq_ = C_seq * E_seq_;
    }
  }

public:
  //! State-space model list
  ListType<std::shared_ptr<_StateSpaceModel>> model_list_;

  //! Total state dimension
  int total_state_dim_ = 0;

  //! Total input dimension
  int total_input_dim_ = 0;

  //! Total output dimension
  int total_output_dim_ = 0;

  //! Sequential extension matrix \f$\boldsymbol{\hat{A}}\f$
  Eigen::Matrix<double, Eigen::Dynamic, StateDim> A_seq_;

  //! Sequential extension matrix \f$\boldsymbol{\hat{B}}\f$
  Eigen::MatrixXd B_seq_;

  //! Sequential extension vector (i.e., offset vector) \f$\boldsymbol{\hat{e}}\f$
  Eigen::VectorXd E_seq_;
};
} // namespace CCC
