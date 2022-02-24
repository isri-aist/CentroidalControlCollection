/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <iostream>

#include <CCC/EigenTypes.h>
#include <CCC/StateSpaceModel.h>

class ModelFixed1 : public CCC::StateSpaceModel<3, 1, 1>
{
public:
  ModelFixed1()
  {
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_(2) = 1;

    C_(0) = 1;
  }
};

class ModelFixed2 : public ModelFixed1
{
public:
  ModelFixed2()
  {
    D_ << 5.0;

    E_ << -1.0, 2.0, -3.0;
  }
};

class ModelDynamic1 : public CCC::StateSpaceModel<3, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  ModelDynamic1() : StateSpaceModel(3, 1, 1)
  {
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_(2) = 1;

    C_(0) = 1;
  }
};

class ModelDynamic2 : public ModelDynamic1
{
public:
  ModelDynamic2()
  {
    D_ << 5.0;

    E_ << -1.0, 2.0, -3.0;
  }
};

template<class ModelType>
void testStateSpaceModel(bool debug = false)
{
  ModelType model;

  Eigen::Vector3d x;
  x << 1, 2, 3;
  Eigen::Vector1d u;
  u << -0.5;

  // Continuous system
  Eigen::Vector3d dx = model.stateEq(x, u);
  Eigen::Vector1d y = model.observEq(x);

  // Discrete system
  double dt = 1e-2;
  model.calcDiscMatrix(dt);
  Eigen::Vector3d nextXDisc = model.stateEqDisc(x, u);

  // Check results
  if(debug)
  {
    std::cout << "x: " << x.transpose() << std::endl;
    std::cout << "dx: " << dx.transpose() << std::endl;
    std::cout << "y: " << y.transpose() << std::endl;

    std::cout << "A:\n" << model.A_ << std::endl;
    std::cout << "B:\n" << model.B_ << std::endl;
    std::cout << "C:\n" << model.C_ << std::endl;
    if(model.D_.norm() > 0)
    {
      std::cout << "D:\n" << model.D_ << std::endl;
    }
    if(model.E_.norm() > 0)
    {
      std::cout << "E: " << model.E_.transpose() << std::endl;
    }

    std::cout << "dt: " << model.dt_ << std::endl;
    std::cout << "Ad:\n" << model.Ad_ << std::endl;
    std::cout << "Bd:\n" << model.Bd_ << std::endl;
    if(model.E_.norm() > 0)
    {
      std::cout << "Ed: " << model.Ed_.transpose() << std::endl;
    }
  }

  Eigen::Vector3d nextXCont = x + dt * dx;
  if(debug)
  {
    std::cout << "continuous result:\n" << nextXCont.transpose() << std::endl;
    std::cout << "discrete result:\n" << nextXDisc.transpose() << std::endl;
  }

  EXPECT_TRUE(nextXCont.isApprox(nextXDisc, 1e-3));
}

TEST(TestStateSpaceModel, Fixed1)
{
  testStateSpaceModel<ModelFixed1>();
}
TEST(TestStateSpaceModel, Fixed2)
{
  testStateSpaceModel<ModelFixed2>();
}
TEST(TestStateSpaceModel, Dynamic1)
{
  testStateSpaceModel<ModelDynamic1>();
}
TEST(TestStateSpaceModel, Dynamic2)
{
  testStateSpaceModel<ModelDynamic2>();
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
