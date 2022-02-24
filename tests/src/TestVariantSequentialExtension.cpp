/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <iostream>

#include <CCC/EigenTypes.h>
#include <CCC/VariantSequentialExtension.h>

class Model1 : public CCC::StateSpaceModel<3, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  Model1() : StateSpaceModel(3, 1, 1)
  {
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_(2) = 1;

    C_(0) = 1;
  }
};

class Model2 : public Model1
{
public:
  Model2()
  {
    D_ << 5.0;

    E_ << -1.0, 2.0, -3.0;
  }
};

class Model3 : public CCC::StateSpaceModel<3, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  Model3(int idx) : StateSpaceModel(3, 0, 1)
  {
    A_(0, 1) = 0.1 * idx;
    A_(1, 2) = 0.2 * idx * idx;

    C_(0) = 1;

    E_ << -0.1 * idx, 0.2 * idx * idx, -0.3;
  }
};

template<class ModelType>
void testVariantSequentialExtension(bool debug = false)
{
  constexpr int state_dim = 3;
  size_t seq_len = 5;

  // Setup model
  std::vector<std::shared_ptr<CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>>> model_list;
  model_list.assign(seq_len, std::make_shared<ModelType>());
  double dt = 1e-2;
  for(const auto & model : model_list)
  {
    model->calcDiscMatrix(dt);
  }

  // Calculate sequential extension
  CCC::VariantSequentialExtension<state_dim> seq_ext_model(model_list);
  if(debug)
  {
    std::cout << "A_seq:\n" << seq_ext_model.A_seq_ << std::endl;
    std::cout << "B_seq:\n" << seq_ext_model.B_seq_ << std::endl;
    std::cout << "E_seq:\n" << seq_ext_model.E_seq_ << std::endl;
  }

  Eigen::Vector3d x_initial;
  x_initial << 1, 2, 3;
  Eigen::VectorXd u_seq(seq_len);
  u_seq << 5.0, 2.5, 0.0, -1.0, -2.0;

  Eigen::VectorXd x_seq_ext = seq_ext_model.A_seq_ * x_initial + seq_ext_model.B_seq_ * u_seq + seq_ext_model.E_seq_;

  // Compare results from extension and iteration
  Eigen::VectorXd x_seq_iter(state_dim * seq_len);
  for(size_t i = 0; i < seq_len; i++)
  {
    if(i == 0)
    {
      x_seq_iter.template segment<state_dim>(i * state_dim) =
          model_list[i]->stateEqDisc(x_initial, u_seq.segment(i, 1));
    }
    else
    {
      x_seq_iter.template segment<state_dim>(i * state_dim) =
          model_list[i]->stateEqDisc(x_seq_iter.template segment<state_dim>((i - 1) * state_dim), u_seq.segment(i, 1));
    }
  }

  if(debug)
  {
    Eigen::Map<Eigen::MatrixXd> x_seq_ext_reshaped(x_seq_ext.data(), state_dim, seq_len);
    Eigen::Map<Eigen::MatrixXd> x_seq_iter_reshaped(x_seq_iter.data(), state_dim, seq_len);
    std::cout << "x_seq_ext:\n" << x_seq_ext_reshaped << std::endl;
    std::cout << "x_seq_iter:\n" << x_seq_iter_reshaped << std::endl;
    std::cout << "error:\n" << x_seq_ext_reshaped - x_seq_iter_reshaped << std::endl;
  }

  EXPECT_TRUE((x_seq_ext - x_seq_iter).norm() < 1e-10);
}

TEST(TestVariantSequentialExtension, Model1)
{
  testVariantSequentialExtension<Model1>();
}

TEST(TestVariantSequentialExtension, Model2)
{
  testVariantSequentialExtension<Model2>();
}

TEST(TestVariantSequentialExtension, VariantInput)
{
  bool debug = false;
  constexpr int state_dim = 3;
  size_t seq_len = 10;
  std::set<size_t> no_input_steps = {4, 5, 6, 8};

  // Setup model
  std::vector<std::shared_ptr<CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>>> model_list(seq_len);
  for(size_t i = 0; i < seq_len; i++)
  {
    if(no_input_steps.count(i))
    {
      model_list[i] = std::make_shared<Model3>(i);
    }
    else
    {
      model_list[i] = std::make_shared<Model1>();
    }
  }
  double dt = 1e-2;
  for(const auto & model : model_list)
  {
    model->calcDiscMatrix(dt);
  }

  // Calculate sequential extension
  CCC::VariantSequentialExtension<state_dim> seq_ext_model(model_list);
  if(debug)
  {
    std::cout << "A_seq:\n" << seq_ext_model.A_seq_ << std::endl;
    std::cout << "B_seq:\n" << seq_ext_model.B_seq_ << std::endl;
    std::cout << "E_seq:\n" << seq_ext_model.E_seq_ << std::endl;
  }

  Eigen::Vector3d x_initial;
  x_initial << 1, 2, 3;
  Eigen::VectorXd u_seq(seq_len - no_input_steps.size());
  u_seq << 0.5, 1.0, 0.5, 1.0, -1.0, -2.0;

  Eigen::VectorXd x_seq_ext = seq_ext_model.A_seq_ * x_initial + seq_ext_model.B_seq_ * u_seq + seq_ext_model.E_seq_;

  // Compare results from extension and iteration
  Eigen::VectorXd x_seq_iter(state_dim * seq_len);
  size_t accum_input_dim = 0;
  for(size_t i = 0; i < seq_len; i++)
  {
    int current_input_dim = (no_input_steps.count(i) ? 0 : 1);
    if(i == 0)
    {
      x_seq_iter.template segment<state_dim>(i * state_dim) =
          model_list[i]->stateEqDisc(x_initial, u_seq.segment(accum_input_dim, current_input_dim));
    }
    else
    {
      x_seq_iter.template segment<state_dim>(i * state_dim) =
          model_list[i]->stateEqDisc(x_seq_iter.template segment<state_dim>((i - 1) * state_dim),
                                     u_seq.segment(accum_input_dim, current_input_dim));
    }
    accum_input_dim += current_input_dim;
  }

  if(debug)
  {
    Eigen::Map<Eigen::MatrixXd> x_seq_ext_reshaped(x_seq_ext.data(), state_dim, seq_len);
    Eigen::Map<Eigen::MatrixXd> x_seq_iter_reshaped(x_seq_iter.data(), state_dim, seq_len);
    std::cout << "x_seq_ext:\n" << x_seq_ext_reshaped << std::endl;
    std::cout << "x_seq_iter:\n" << x_seq_iter_reshaped << std::endl;
    std::cout << "error:\n" << x_seq_ext_reshaped - x_seq_iter_reshaped << std::endl;
  }

  EXPECT_TRUE((x_seq_ext - x_seq_iter).norm() < 1e-10);
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
