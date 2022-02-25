/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <fstream>

#include <CCC/LinearMpcTZ.h>

TEST(TestLinearMpcTZ, Test1)
{
  double mass = 100.0; // [kg]
  double dt = 0.05; // [ms]
  double motion_duration = 10.0; // [sec]

  int seq_len = static_cast<int>(motion_duration / dt);
  std::vector<bool> contact_seq(seq_len);
  Eigen::VectorXd ref_pos_seq(seq_len);
  for(int i = 0; i < seq_len; i++)
  {
    double t = i * dt;

    if(5.0 < t && t < 5.25 || 6.0 < t && t < 6.5)
    {
      contact_seq[i] = false;
    }
    else
    {
      contact_seq[i] = true;
    }

    if(t < 8.5)
    {
      ref_pos_seq[i] = 1.0;
    }
    else
    {
      ref_pos_seq[i] = 0.8;
    }
  }
  double initial_pos = 0.8; // [m]
  double initial_vel = -1.0; // [m/s]

  CCC::LinearMpcTZ mpc(mass, dt);
  mpc.runLoop(contact_seq, initial_pos, initial_vel, ref_pos_seq);

  std::ofstream ofs("/tmp/TestLinearMpcTZ.txt");
  ofs << "time contact ref_pos planned_pos planned_vel planned_force" << std::endl;
  for(int i = 0; i < seq_len; i++)
  {
    ofs << i * dt << " " << contact_seq[i] << " " << ref_pos_seq[i] << " " << mpc.planned_pos_seq_[i] << " "
        << mpc.planned_vel_seq_[i] << " " << mpc.planned_force_seq_[i] << std::endl;
  }
  std::cout << "Run the following commands in gnuplot:\n"
               "  set key autotitle columnhead\n"
               "  set key noenhanced\n"
               "  plot \"/tmp/TestLinearMpcTZ.txt\" u 1:3 w lp, \"\" u 1:4 w lp\n";

  EXPECT_LT(std::abs(mpc.planned_pos_seq_[seq_len - 1] - ref_pos_seq[seq_len - 1]), 1e-2);
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
