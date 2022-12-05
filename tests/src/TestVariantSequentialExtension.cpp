/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <deque>
#include <iostream>

#include <CCC/VariantSequentialExtension.h>

constexpr int state_dim = 3;

class Model1 : public CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  Model1() : StateSpaceModel(state_dim, 1, 1)
  {
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_(2) = 1;

    C_(0, 0) = 1;
    C_(0, 2) = 2;
  }
};

class Model2 : public Model1
{
public:
  Model2()
  {
    E_ << -1.0, 2.0, -3.0;
  }
};

class Model3 : public CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  Model3(int idx) : StateSpaceModel(state_dim, 0, 1)
  {
    A_(0, 1) = 0.1 * idx;
    A_(1, 2) = 0.2 * idx * idx;

    C_(0, 0) = 1;
    C_(0, 1) = -1;

    E_ << -0.1 * idx, 0.2 * idx * idx, -0.3;
  }
};

template<template<class, class...> class ListType>
Eigen::VectorXd calcStateSeq(const CCC::VariantSequentialExtension<state_dim, ListType> & seq_ext,
                             const Eigen::VectorXd & u_seq,
                             const Eigen::Matrix<double, state_dim, 1> & x_initial,
                             bool extend_for_output)
{
  Eigen::VectorXd x_seq_iter(extend_for_output ? seq_ext.totalOutputDim() : seq_ext.totalStateDim());
  size_t accum_input_dim = 0;
  size_t accum_output_dim = 0;
  Eigen::Matrix<double, state_dim, 1> current_x = x_initial;
  for(size_t i = 0; i < seq_ext.model_list_.size(); i++)
  {
    int current_input_dim = seq_ext.model_list_[i]->inputDim();
    int current_output_dim = seq_ext.model_list_[i]->outputDim();
    current_x = seq_ext.model_list_[i]->stateEqDisc(current_x, u_seq.segment(accum_input_dim, current_input_dim));
    if(extend_for_output)
    {
      x_seq_iter.segment(accum_output_dim, current_output_dim) = seq_ext.model_list_[i]->observEq(current_x);
    }
    else
    {
      x_seq_iter.template segment<state_dim>(i * state_dim) = current_x;
    }
    accum_input_dim += current_input_dim;
    accum_output_dim += current_output_dim;
  }
  return x_seq_iter;
}

template<template<class, class...> class ListType>
void printStateSeq(const CCC::VariantSequentialExtension<state_dim, ListType> & seq_ext,
                   Eigen::VectorXd x_seq_ext, // Pass by value for Eigen::Map
                   Eigen::VectorXd x_seq_iter,
                   bool extend_for_output)
{
  if(extend_for_output)
  {
    {
      std::cout << "x_seq_ext:\n";
      int accum_output_dim = 0;
      for(size_t i = 0; i < seq_ext.model_list_.size(); i++)
      {
        int current_output_dim = seq_ext.model_list_[i]->outputDim();
        std::cout << x_seq_ext.segment(accum_output_dim, current_output_dim) << std::endl;
        accum_output_dim += current_output_dim;
      }
    }
    {
      std::cout << "x_seq_iter:\n";
      int accum_output_dim = 0;
      for(size_t i = 0; i < seq_ext.model_list_.size(); i++)
      {
        int current_output_dim = seq_ext.model_list_[i]->outputDim();
        std::cout << x_seq_iter.segment(accum_output_dim, current_output_dim) << std::endl;
        accum_output_dim += current_output_dim;
      }
    }
    {
      std::cout << "error:\n";
      int accum_output_dim = 0;
      for(size_t i = 0; i < seq_ext.model_list_.size(); i++)
      {
        int current_output_dim = seq_ext.model_list_[i]->outputDim();
        std::cout << (x_seq_ext - x_seq_iter).segment(accum_output_dim, current_output_dim) << std::endl;
        accum_output_dim += current_output_dim;
      }
    }
  }
  else
  {
    Eigen::Map<Eigen::MatrixXd> x_seq_ext_reshaped(x_seq_ext.data(), state_dim, seq_ext.model_list_.size());
    Eigen::Map<Eigen::MatrixXd> x_seq_iter_reshaped(x_seq_iter.data(), state_dim, seq_ext.model_list_.size());
    std::cout << "x_seq_ext:\n" << x_seq_ext_reshaped.transpose() << std::endl;
    std::cout << "x_seq_iter:\n" << x_seq_iter_reshaped.transpose() << std::endl;
    std::cout << "error:\n" << (x_seq_ext_reshaped - x_seq_iter_reshaped).transpose() << std::endl;
  }
}

template<class ModelType>
void testVariantSequentialExtension(bool extend_for_output, bool debug = false)
{
  // Setup model
  size_t seq_len = 5;
  double dt = 1e-2;
  const auto & model = std::make_shared<ModelType>();
  model->calcDiscMatrix(dt);
  std::vector<std::shared_ptr<CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>>> model_list(seq_len,
                                                                                                           model);

  // Calculate sequential extension
  CCC::VariantSequentialExtension<state_dim> seq_ext(model_list, extend_for_output);
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

TEST(TestVariantSequentialExtension, Model1)
{
  testVariantSequentialExtension<Model1>(false);
  testVariantSequentialExtension<Model1>(true);
}

TEST(TestVariantSequentialExtension, Model2)
{
  testVariantSequentialExtension<Model2>(false);
  testVariantSequentialExtension<Model2>(true);
}

void testVariantSequentialExtensionVariantInput(bool extend_for_output, bool debug = false)
{
  // Setup model
  size_t seq_len = 10;
  double dt = 1e-2;
  std::set<size_t> no_input_steps = {4, 5, 6, 8};
  std::deque<std::shared_ptr<CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>>> model_list;
  for(size_t i = 0; i < seq_len; i++)
  {
    std::shared_ptr<CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>> model;
    if(no_input_steps.count(i))
    {
      model = std::make_shared<Model3>(i);
    }
    else
    {
      model = std::make_shared<Model1>();
    }
    model->calcDiscMatrix(dt);
    model_list.push_back(model);
  }

  // Calculate sequential extension
  CCC::VariantSequentialExtension<state_dim, std::deque> seq_ext(model_list, extend_for_output);
  if(debug)
  {
    std::cout << "A_seq:\n" << seq_ext.A_seq_ << std::endl;
    std::cout << "B_seq:\n" << seq_ext.B_seq_ << std::endl;
    std::cout << "E_seq:\n" << seq_ext.E_seq_ << std::endl;
  }

  // Calculate state sequence
  Eigen::Matrix<double, state_dim, 1> x_initial;
  x_initial << 1, 2, 3;
  Eigen::VectorXd u_seq(seq_len - no_input_steps.size());
  u_seq << 0.5, 1.0, 0.5, 1.0, -1.0, -2.0;

  Eigen::VectorXd x_seq_ext = seq_ext.A_seq_ * x_initial + seq_ext.B_seq_ * u_seq + seq_ext.E_seq_;
  Eigen::VectorXd x_seq_iter = calcStateSeq(seq_ext, u_seq, x_initial, extend_for_output);

  // Compare results from extension and iteration
  EXPECT_LT((x_seq_ext - x_seq_iter).norm(), 1e-10);
  if(debug)
  {
    printStateSeq(seq_ext, x_seq_ext, x_seq_iter, extend_for_output);
  }
}

TEST(TestVariantSequentialExtension, VariantInput)
{
  testVariantSequentialExtensionVariantInput(false);
  testVariantSequentialExtensionVariantInput(true);
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
