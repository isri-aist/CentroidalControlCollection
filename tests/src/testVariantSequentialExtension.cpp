/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <iostream>

#include <CCC/EigenTypes.h>
#include <CCC/VariantSequentialExtension.h>

class TestModel1 : public CCC::StateSpaceModel<3, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  TestModel1() : StateSpaceModel(3, 1, 1)
  {
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_(2) = 1;

    C_(0) = 1;
  }
};

class TestModel2 : public TestModel1
{
public:
  TestModel2()
  {
    D_ << 5.0;

    E_ << -1.0, 2.0, -3.0;
  }
};

template<class ModelType>
void testVariantSequentialExtension(bool debug = false)
{
  constexpr int state_dim = 3;
  std::vector<std::shared_ptr<CCC::StateSpaceModel<state_dim, Eigen::Dynamic, Eigen::Dynamic>>> model_list;
  model_list.assign(5, std::make_shared<ModelType>());

  double dt = 1e-2;
  for(const auto & model : model_list)
  {
    model->calcDiscMatrix(dt);
  }

  CCC::VariantSequentialExtension<state_dim> seq_ext_model(model_list);

  std::cout << "A_seq:\n" << seq_ext_model.A_seq_ << std::endl;
  std::cout << "B_seq:\n" << seq_ext_model.B_seq_ << std::endl;
  std::cout << "E_seq:\n" << seq_ext_model.E_seq_ << std::endl;
}

TEST(TestVariantSequentialExtension, Test1)
{
  testVariantSequentialExtension<TestModel1>();
}

TEST(TestVariantSequentialExtension, Test2)
{
  testVariantSequentialExtension<TestModel2>();
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
