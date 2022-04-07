/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include <CCC/Constants.h>
#include <CCC/DdpZmp.h>
#include <CCC/EigenTypes.h>
#include <CCC/StateSpaceModel.h>

/** \brief State-space model of CoM-ZMP dynamics with ZMP input. */
class ComZmpSimModel : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
  */
  ComZmpSimModel(double com_height)
  {
    double omega2 = CCC::constants::g / com_height;

    A_(0, 1) = 1;
    A_(1, 0) = omega2;

    B_(1, 0) = -1 * omega2;
  }
};

struct ComPosVel
{
  Eigen::Vector2d x = Eigen::Vector2d::Zero();
  Eigen::Vector2d y = Eigen::Vector2d::Zero();
  Eigen::Vector2d z = Eigen::Vector2d::Zero();
};

/** \brief State-space model of vertical motion dynamics with force input. */
class VerticalSimModel : public CCC::StateSpaceModel<2, 1, 0>
{
public:
  /** \brief Constructor.
      \param mass robot mass [Kg]
  */
  VerticalSimModel(double mass)
  {
    A_(0, 1) = 1;

    B_(1, 0) = 1 / mass;

    E_(1) = -1 * CCC::constants::g;
  }
};

TEST(TestDdpZmp, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.005; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double mass = 100.0; // [Kg]
  double ref_com_height = 1.0; // [m]

  // Setup DDP
  CCC::DdpZmp ddp(mass, horizon_dt, horizon_steps);

  std::function<double(double)> ref_zmp_func = [](double t) {
    if(t < 2.0)
    {
      return 0.0;
    }
    else if(t < 3.0)
    {
      return 0.1;
    }
    else if(t < 4.0)
    {
      return -0.1;
    }
    else if(t < 5.0)
    {
      return 0.1;
    }
    else if(t < 6.0)
    {
      return 0.2;
    }
    else if(t < 7.0)
    {
      return 0.3;
    }
    else
    {
      return 0.4;
    }
  };
  std::function<CCC::DdpZmp::RefData(double)> ref_data_func = [&](double t) {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    CCC::DdpZmp::RefData ref_data;
    ref_data.zmp.setConstant(ref_zmp_func(t));
    ref_data.com_z = ref_com_height;

    return ref_data;
  };

  // Run control
  ComPosVel com_pos_vel;
  com_pos_vel.z[0] = ref_com_height;
  CCC::DdpZmp::PlannedData planned_data;

  std::string file_path = "/tmp/TestDdpZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos com_vel planned_zmp ref_zmp" << std::endl;

  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    CCC::DdpZmp::InitialParam initial_param;
    initial_param.pos << com_pos_vel.x[0], com_pos_vel.y[0], com_pos_vel.z[0];
    initial_param.vel << com_pos_vel.x[1], com_pos_vel.y[1], com_pos_vel.z[1];
    planned_data = ddp.planOnce(ref_data_func, initial_param, t);

    // Dump
    const auto & ref_data = ref_data_func(t);
    ofs << t << " " << com_pos_vel.x.transpose() << " " << planned_data.zmp.x() << " " << ref_data.zmp.x() << std::endl;

    // Check
    EXPECT_LT(std::abs(planned_data.zmp.x() - ref_data.zmp.x()), 0.1); // [m]

    // Setup simulation
    ComZmpSimModel sim_model_xy(com_pos_vel.z[0]);
    sim_model_xy.calcDiscMatrix(horizon_dt);
    VerticalSimModel sim_model_z(mass);
    sim_model_z.calcDiscMatrix(horizon_dt);

    // Update
    t += horizon_dt;
    Eigen::Vector1d sim_input;
    sim_input << planned_data.zmp.x();
    com_pos_vel.x = sim_model_xy.stateEqDisc(com_pos_vel.x, sim_input);
    sim_input << planned_data.zmp.y();
    com_pos_vel.y = sim_model_xy.stateEqDisc(com_pos_vel.y, sim_input);
    sim_input << planned_data.force_z;
    com_pos_vel.z = sim_model_z.stateEqDisc(com_pos_vel.z, sim_input);
  }

  // Final check
  const auto & ref_data = ref_data_func(t);
  EXPECT_LT(std::abs(planned_data.zmp.x() - ref_data.zmp.x()), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel.x[0] - ref_data.zmp.x()), 1e-2);
  EXPECT_LT(std::abs(com_pos_vel.x[1]), 1e-2);

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp, \"\" u 1:5 w l lw 2\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
