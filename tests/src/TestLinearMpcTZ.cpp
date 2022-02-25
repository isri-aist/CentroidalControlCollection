/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <CCC/LinearMpcTZ.h>

TEST(TestLinearMpcTZ, Test1)
{
  double mass = 100.0; // [kg]
  double horizon_dt = 0.05; // [ms]

  CCC::LinearMpcTZ mpc(mass, horizon_dt);

  std::function<bool(double)> contact_func = [](double t) { return !(5.0 < t && t < 5.25 || 6.0 < t && t < 6.5); };
  std::function<double(double)> ref_pos_func = [](double t) { return t < 8.5 ? 1.0 : 0.8; }; // [m]
  Eigen::Vector2d initial_pos_vel(0.8, -1.0); // ([m], [m/s])
  std::pair<double, double> motion_time_range(0.0, 10.0); // ([s], [s])
  double horizon_duration = 4.0; // [sec]
  double sim_dt = 0.04; // [ms]

  mpc.planLoop(contact_func, ref_pos_func, initial_pos_vel, motion_time_range, horizon_duration, sim_dt);

  mpc.dumpResultDataSeq("/tmp/TestLinearMpcTZ.txt", true);

  EXPECT_LT(std::abs(mpc.result_data_seq_[mpc.result_data_seq_.size() - 1].planned_pos
                     - ref_pos_func(motion_time_range.second)),
            1e-2); // [m]
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
