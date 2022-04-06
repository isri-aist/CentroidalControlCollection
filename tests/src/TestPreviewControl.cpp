/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/PreviewControl.h>

class TestBipedalModel : public CCC::StateSpaceModel<3, 1, 1>
{
public:
  TestBipedalModel(double comHeight = 1.0)
  {
    A_ = Eigen::Matrix<double, 3, 3>::Zero();
    A_(0, 1) = 1;
    A_(1, 2) = 1;

    B_ = Eigen::Matrix<double, 3, 1>::Zero();
    B_(2, 0) = 0.1;

    C_ = Eigen::Matrix<double, 1, 3>::Zero();
    C_(0, 0) = 1.0;
    C_(0, 2) = -comHeight / CCC::constants::g;
  }
};

TEST(TestPreviewControl, Test1)
{
  std::shared_ptr<TestBipedalModel> model = std::make_shared<TestBipedalModel>();

  CCC::PreviewControl<3, 1, 1> pc(model);

  // Calculate gain
  double horizon = 2.0; // [sec]
  double dt = 0.005; // [sec]
  Eigen::Vector1d outputWeight;
  outputWeight << 1.0;
  Eigen::Vector1d inputWeight;
  inputWeight << 1e-8;
  pc.calcGain(horizon, dt, outputWeight, inputWeight);

  // Run control
  Eigen::Vector3d state = Eigen::Vector3d::Zero();
  Eigen::Vector1d input;
  Eigen::Vector1d output;
  Eigen::VectorXd refOutputTraj(model->outputDim() * pc.previewSize_);

  std::string file_path = "/tmp/TestPreviewControl.txt";
  std::ofstream ofs(file_path);
  ofs << "# Time, State(" << model->stateDim() << "), Input(" << model->inputDim() << "), RefOutput("
      << model->outputDim() << "), PlannedOutput(" << model->outputDim() << ")" << std::endl;

  double t = 0;
  while(t < 10.0)
  {
    // Generate reference trajectory
    for(int i = 0; i < pc.previewSize_; i++)
    {
      double _t = t + i * dt;
      if(_t < 2.0)
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()).setZero();
      }
      else if(_t < 3.0)
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()) << 0.1;
      }
      else if(_t < 4.0)
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()) << -0.1;
      }
      else if(_t < 5.0)
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()) << 0.1;
      }
      else if(_t < 6.0)
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()) << 0.2;
      }
      else if(_t < 7.0)
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()) << 0.3;
      }
      else
      {
        refOutputTraj.segment(i * model->outputDim(), model->outputDim()) << 0.4;
      }
    }

    // Calculate input
    input = pc.calcOptimalInput(state, refOutputTraj);
    output = model->observEq(state);

    // Dump
    ofs << t << " " << state.transpose() << " " << input.transpose() << " "
        << refOutputTraj.head(model->outputDim()).transpose() << " " << output.transpose() << std::endl;

    // Update time and state
    t += dt;
    state = model->stateEqDisc(state, input);
  }

  std::cout << "Plot the results by the following commands in gnuplot." << std::endl;
  std::cout << "  plot \"" << file_path
            << "\" u 1:2 w lp t \"CoM\", \"\" u 1:6 w lp t \"ref ZMP\", \"\" u 1:7 w lp t \"planned ZMP\"" << std::endl;
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
