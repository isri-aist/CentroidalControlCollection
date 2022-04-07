/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/PreviewControlZmp.h>

/** \brief State-space model of CoM-ZMP dynamics with ZMP input. */
class ComZmpSimModel : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
  */
  ComZmpSimModel(double com_height)
  {
    A_(0, 1) = 1;
    A_(1, 0) = CCC::constants::g / com_height;

    B_(1, 0) = -1 * CCC::constants::g / com_height;
  }
};

TEST(TestPreviewControlZmp, Test1)
{
  double com_height = 1.0; // [m]
  CCC::PreviewControlZmp pc(std::make_shared<CCC::PreviewControlZmp::ComZmpModel>(com_height));

  // Calculate gain
  double horizon = 2.0; // [sec]
  double dt = 0.005; // [sec]
  pc.calcGain(horizon, dt);

  // Setup simulation model
  ComZmpSimModel sim_model(com_height);
  sim_model.calcDiscMatrix(dt);

  // Run control
  Eigen::Vector2d state = Eigen::Vector2d::Zero();
  double planned_zmp = state[0];
  Eigen::VectorXd ref_zmp_traj(pc.preview_size_);

  std::string file_path = "/tmp/TestPreviewControlZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos com_vel com_acc planned_zmp ref_zmp" << std::endl;

  double t = 0;
  while(t < 10.0)
  {
    // Generate reference trajectory
    for(int i = 0; i < pc.preview_size_; i++)
    {
      double _t = t + i * dt;
      if(_t < 2.0)
      {
        ref_zmp_traj[i] = 0.0;
      }
      else if(_t < 3.0)
      {
        ref_zmp_traj[i] = 0.1;
      }
      else if(_t < 4.0)
      {
        ref_zmp_traj[i] = -0.1;
      }
      else if(_t < 5.0)
      {
        ref_zmp_traj[i] = 0.1;
      }
      else if(_t < 6.0)
      {
        ref_zmp_traj[i] = 0.2;
      }
      else if(_t < 7.0)
      {
        ref_zmp_traj[i] = 0.3;
      }
      else
      {
        ref_zmp_traj[i] = 0.4;
      }
    }

    // Calculate input
    CCC::PreviewControlZmp::InitialParam initial_param;
    initial_param << state, CCC::constants::g / com_height * (state[0] - planned_zmp);
    planned_zmp = pc.planZmp(initial_param, ref_zmp_traj);

    // Dump
    ofs << t << " " << initial_param.transpose() << " " << planned_zmp << " " << ref_zmp_traj[0] << std::endl;

    // Update time and state
    t += dt;
    Eigen::Vector1d input;
    input << planned_zmp;
    state = sim_model.stateEqDisc(state, input);
  }

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:5 w lp, \"\" u 1:6 w l lw 2\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
