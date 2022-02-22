/* Author: Masaki Murooka */

#pragma once

#include <memory>

#include <CCC/StateSpaceModel.h>

namespace CCC
{
/*! \brief Sequential extended model for time-variant system.
 *
 *  \tparam StateDim state dimension
 *
 *  \note State dimension must be the same for all models in the sequence.
 */
template<int StateDim>
class SequentialExtendedModelVariant
{
public:
  /*! \brief Type of state-space model with fixed state dimension. */
  using _StateSpaceModel = StateSpaceModel<StateDim, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*! \brief Constructor.
   *
   *  \param model_list list of state-space model
   */
  SequentialExtendedModelVariant(const std::vector<std::shared_ptr<_StateSpaceModel>>& model_list):
      model_list_(model_list)
  {
    setup();
  }

 protected:
  /*! \brief Setup coefficients. */
  void setup()
  {
    size_t seq_len = model_list_.size();
    int total_input_dim = 0;
    for (const auto& model : model_list_) {
      total_input_dim += model->inputDim();
    }

    // Resize matrix
    A_seq_.setZero(seq_len * StateDim, StateDim);
    B_seq_.setZero(seq_len * StateDim, total_input_dim);
    E_seq_.setZero(seq_len * StateDim);

    int current_total_input_dim = 0;
    for (size_t i = 0; i < seq_len; i++) {
      int current_input_dim = model_list_[i]->inputDim();

      // Setq A_seq_
      if (i == 0) {
        A_seq_.template middleRows<StateDim>(i * StateDim) = model_list_[0]->Ad_;
      } else {
        (A_seq_.template middleRows<StateDim>(i * StateDim)).noalias() =
            model_list_[i]->Ad_ * A_seq_.template middleRows<StateDim>((i - 1) * StateDim);
      }

      // Setq B_seq_
      for (size_t j = i; j < seq_len; j++) {
        if (j == i) {
          B_seq_.block(j * StateDim, current_total_input_dim, StateDim, current_input_dim) = model_list_[i]->Bd_;
        } else {
          (B_seq_.block(j * StateDim, current_total_input_dim, StateDim, current_input_dim)).noalias() =
              model_list_[j]->Ad_ * B_seq_.block((j - 1) * StateDim, current_total_input_dim, StateDim, current_input_dim);
        }
      }

      // Setq E_seq_
      if (i == 0) {
        E_seq_.segment<StateDim>(i * StateDim) = model_list_[0]->Ed_;
      } else {
        E_seq_.segment<StateDim>(i * StateDim) =
            model_list_[i]->Ad_ * E_seq_.segment<StateDim>((i - 1) * StateDim) + model_list_[i]->Ed_;
      }

      current_total_input_dim += model_list_[i]->inputDim();
    }
  }

 public:
  //! State-space model
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list_;

  //! Sequential extended matrix of A in discrete system
  Eigen::Matrix<double, Eigen::Dynamic, StateDim> A_seq_;

  //! Sequential extended matrix of B in discrete system
  Eigen::MatrixXd B_seq_;

  //! Sequential extended vector of E (i.e., offset vector) in discrete system
  Eigen::VectorXd E_seq_;
};
}
