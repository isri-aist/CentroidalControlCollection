/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <CCC/Constants.h>
#include <CCC/LinearMpcZ.h>

TEST(TestLinearMpcZ, Test1)
{
  double mass = 100.0; // [kg]
  double horizon_dt = 0.05; // [ms]

  CCC::LinearMpcZ mpc(mass, horizon_dt);

  std::function<bool(double)> contact_func = [](double t) { return !(5.0 < t && t < 5.25 || 6.0 < t && t < 6.5); };
  std::function<double(double)> ref_pos_func = [](double t) { return t < 8.5 ? 1.0 : 0.8; }; // [m]
  Eigen::Vector2d initial_pos_vel(0.8, -1.0); // ([m], [m/s])
  std::pair<double, double> motion_time_range(0.0, 10.0); // ([s], [s])
  double horizon_duration = 4.0; // [sec]
  double sim_dt = 0.04; // [ms]

  mpc.planLoop(contact_func, ref_pos_func, initial_pos_vel, motion_time_range, horizon_duration, sim_dt);

  mpc.dumpMotionDataSeq("/tmp/TestLinearMpcZ.txt", true);

  // Check final position
  EXPECT_LT(std::abs(mpc.motion_data_seq_[mpc.motion_data_seq_.size() - 1].planned_pos
                     - ref_pos_func(motion_time_range.second)),
            1e-2); // [m]

  // Check acceleration in the air
  for(const auto & motion_data : mpc.motion_data_seq_)
  {
    if(!contact_func(motion_data.time))
    {
      EXPECT_LT(std::abs(motion_data.planned_acc - -1 * CCC::constants::g), 1e-8); // [m/s^2]
    }
  }
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
