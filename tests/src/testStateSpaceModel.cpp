/* Author: Masaki Murooka */

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <EigenTypes.h>
#include <StateSpaceModel.h>

/*! \brief Test state-space model of fixed size. */
class TestModel1 : public CCC::StateSpaceModel<3, 1, 1>
{
public:
  TestModel1()
  {
    A_.setZero();
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_.setZero();
    B_(2) = 1;

    C_.setZero();
    C_(0) = 1;

    D_.setZero();

    E_.setZero();
  }
};

/*! \brief Test state-space model of fixed size with non-zero D and E coefficients. */
class TestModel2 : public TestModel1
{
public:
  TestModel2()
  {
    D_ << 5.0;

    E_ << -1.0, 2.0, -3.0;
  }
};

/*! \brief Test state-space model of dynamic size. */
class TestModel3 : public CCC::StateSpaceModel<3, Eigen::Dynamic, Eigen::Dynamic>
{
public:
  TestModel3()
  {
    A_.setZero(3, 3);
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_.setZero(3, 1);
    B_(2) = 1;

    C_.setZero(1, 3);
    C_(0) = 1;
  }
};

template <class ModelType>
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
  if (debug) {
    std::cout << "x: " << x.transpose() << std::endl;
    std::cout << "dx: " << dx.transpose() << std::endl;
    std::cout << "y: " << y.transpose() << std::endl;

    std::cout << "A:\n" << model.A_ << std::endl;
    std::cout << "B:\n" << model.B_ << std::endl;
    std::cout << "C:\n" << model.C_ << std::endl;
    if (model.D_.norm() > 0) {
      std::cout << "D:\n" << model.D_ << std::endl;
    }
    if (model.E_.norm() > 0) {
      std::cout << "E: " << model.E_.transpose() << std::endl;
    }

    std::cout << "dt: " << model.dt_ << std::endl;
    std::cout << "Ad:\n" << model.Ad_ << std::endl;
    std::cout << "Bd:\n" << model.Bd_ << std::endl;
    if (model.E_.norm() > 0) {
      std::cout << "Ed: " << model.Ed_.transpose() << std::endl;
    }
  }

  Eigen::Vector3d nextXCont = x + dt * dx;
  if (debug) {
    std::cout << "continuous result:\n" << nextXCont.transpose() << std::endl;
    std::cout << "discrete result:\n" << nextXDisc.transpose() << std::endl;
  }

  BOOST_CHECK(nextXCont.isApprox(nextXDisc, 1e-3));
}

BOOST_AUTO_TEST_CASE(TestStateSpaceModel1) { testStateSpaceModel<TestModel1>(); }
BOOST_AUTO_TEST_CASE(TestStateSpaceModel2) { testStateSpaceModel<TestModel2>(); }
BOOST_AUTO_TEST_CASE(TestStateSpaceModel3) { testStateSpaceModel<TestModel3>(); }
