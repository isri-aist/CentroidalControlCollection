/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <iostream>

#include <CCC/InvariantSequentialExtension.h>
#include <CCC/VariantSequentialExtension.h>

constexpr int state_dim = 3;
constexpr int input_dim = 1;
constexpr int output_dim = 1;

class Model1 : public CCC::StateSpaceModel<state_dim, input_dim, output_dim>
{
public:
  Model1() : StateSpaceModel(state_dim, input_dim, output_dim)
  {
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_(2) = 1;

    C_(0, 0) = 1;
    C_(0, 2) = 2;
  }
};

Eigen::VectorXd calcStateSeq(const CCC::InvariantSequentialExtension<state_dim, input_dim, output_dim> & seq_ext,
                             const Eigen::VectorXd & u_seq,
                             const Eigen::Matrix<double, state_dim, 1> & x_initial,
                             bool extend_for_output)
{
  Eigen::VectorXd x_seq_iter(extend_for_output ? seq_ext.totalOutputDim() : seq_ext.totalStateDim());
  Eigen::Matrix<double, state_dim, 1> current_x = x_initial;
  for(size_t i = 0; i < seq_ext.seq_len_; i++)
  {
    current_x = seq_ext.model_->stateEqDisc(current_x, u_seq.segment<input_dim>(i * input_dim));
    if(extend_for_output)
    {
      x_seq_iter.segment<output_dim>(i * output_dim) = seq_ext.model_->observEq(current_x);
    }
    else
    {
      x_seq_iter.segment<state_dim>(i * state_dim) = current_x;
    }
  }
  return x_seq_iter;
}

void printStateSeq(const CCC::InvariantSequentialExtension<state_dim, input_dim, output_dim> & seq_ext,
                   Eigen::VectorXd x_seq_ext, // Pass by value for Eigen::Map
                   Eigen::VectorXd x_seq_iter,
                   bool extend_for_output)
{
  int row_dim = extend_for_output ? output_dim : state_dim;
  Eigen::Map<Eigen::MatrixXd> x_seq_ext_reshaped(x_seq_ext.data(), row_dim, seq_ext.seq_len_);
  Eigen::Map<Eigen::MatrixXd> x_seq_iter_reshaped(x_seq_iter.data(), row_dim, seq_ext.seq_len_);
  std::cout << "x_seq_ext:\n" << x_seq_ext_reshaped.transpose() << std::endl;
  std::cout << "x_seq_iter:\n" << x_seq_iter_reshaped.transpose() << std::endl;
  std::cout << "error:\n" << (x_seq_ext_reshaped - x_seq_iter_reshaped).transpose() << std::endl;
}

template<class ModelType>
void testInvariantSequentialExtension(bool extend_for_output, bool debug = false)
{
  // Setup model
  size_t seq_len = 5;
  double dt = 1e-2;
  const auto & model = std::make_shared<ModelType>();
  model->calcDiscMatrix(dt);

  // Calculate sequential extension
  CCC::InvariantSequentialExtension<state_dim, input_dim, output_dim> seq_ext(model, seq_len, extend_for_output);
  if(debug)
  {
    std::cout << "A_seq:\n" << seq_ext.A_seq_ << std::endl;
    std::cout << "B_seq:\n" << seq_ext.B_seq_ << std::endl;
    std::cout << "E_seq:\n" << seq_ext.E_seq_ << std::endl;
  }

  // Calculate state sequence
  Eigen::Matrix<double, state_dim, 1> x_initial;
  x_initial << 1, 2, 3;
  Eigen::VectorXd u_seq(seq_len);
  u_seq << 5.0, 2.5, 0.0, -1.0, -2.0;

  Eigen::VectorXd x_seq_ext = seq_ext.A_seq_ * x_initial + seq_ext.B_seq_ * u_seq + seq_ext.E_seq_;
  Eigen::VectorXd x_seq_iter = calcStateSeq(seq_ext, u_seq, x_initial, extend_for_output);

  // Compare results from extension and iteration
  EXPECT_LT((x_seq_ext - x_seq_iter).norm(), 1e-10);
  if(debug)
  {
    printStateSeq(seq_ext, x_seq_ext, x_seq_iter, extend_for_output);
  }
}

TEST(TestInvariantSequentialExtension, Model1)
{
  testInvariantSequentialExtension<Model1>(false);
  testInvariantSequentialExtension<Model1>(true);
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
