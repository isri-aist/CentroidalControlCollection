/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <CCC/Constants.h>
#include <CCC/LinearMpcZ.h>

TEST(TestLinearMpcZ, Test1)
{
  double mass = 100.0; // [kg]
  double horizon_dt = 0.05; // [sec]

  CCC::LinearMpcZ mpc(mass, horizon_dt);

  std::function<bool(double)> contact_func = [](double t) { return !(5.0 < t && t < 5.25 || 6.0 < t && t < 6.5); };
  std::function<double(double)> ref_pos_func = [](double t) { return t < 8.5 ? 1.0 : 0.8; }; // [m]
  CCC::LinearMpcZ::InitialParam initial_param(0.8, -1.0); // ([m], [m/s])
  std::pair<double, double> motion_time_range(0.0, 10.0); // ([sec], [sec])
  double horizon_duration = 2.0; // [sec]
  double sim_dt = 0.04; // [sec]

  mpc.planLoop(contact_func, ref_pos_func, initial_param, motion_time_range, horizon_duration, sim_dt);

  mpc.dumpMotionDataSeq("/tmp/TestLinearMpcZ.txt", true);

  // Check final state
  const auto & motion_data_final = mpc.motionDataSeq().rbegin()->second;
  EXPECT_LT(std::abs(motion_data_final.planned_pos - motion_data_final.ref_pos),
            1e-2); // [m]
  EXPECT_LT(std::abs(motion_data_final.planned_vel),
            1e-2); // [m/s]
  EXPECT_LT(std::abs(motion_data_final.ref_pos - ref_pos_func(motion_time_range.second)),
            1e-10); // [m]

  // Check acceleration in the air
  for(const auto & motion_data_kv : mpc.motionDataSeq())
  {
    if(!contact_func(motion_data_kv.first))
    {
      EXPECT_LT(std::abs(motion_data_kv.second.planned_acc - -1 * CCC::constants::g), 1e-8); // [m/s^2]
    }
  }
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
