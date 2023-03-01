/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <CCC/LinearMpcZ.h>

#include "SimModels.h"

TEST(TestLinearMpcZ, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.05; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.04; // [sec]
  double mass = 100.0; // [Kg]

  // Setup MPC
  std::vector<double> computation_duration_list;
  CCC::LinearMpcZ mpc(mass, horizon_dt, horizon_steps);

  // Setup contact
  std::function<bool(double)> contact_func = [](double t) { return !((5.0 < t && t < 5.25) || (6.0 < t && t < 6.5)); };
  std::function<double(double)> ref_pos_func = [](double t) { return t < 8.5 ? 1.0 : 0.8; }; // [m]

  // Setup simulation
  VerticalSimModel sim(mass, sim_dt);
  sim.state_ << ref_pos_func(0.0), 0.0;

  // Setup dump file
  std::string file_path = "/tmp/TestLinearMpcZ.txt";
  std::ofstream ofs(file_path);
  ofs << "time planned_pos planned_vel ref_pos planned_force computation_time" << std::endl;

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    CCC::LinearMpcZ::InitialParam initial_param = sim.state_;
    double planned_force = mpc.planOnce(contact_func, ref_pos_func, initial_param, t);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    // Dump
    double ref_pos = ref_pos_func(t);
    ofs << t << " " << sim.state_.transpose() << " " << ref_pos << " " << planned_force << " "
        << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT(std::abs(sim.state_[0] - ref_pos), 2.0); // [m]
    EXPECT_LT(std::abs(sim.state_[1]), 5.0); // [m/s]
    if(!contact_func(t))
    {
      EXPECT_LT(std::abs(planned_force), 1e-8); // [N]
    }

    // Simulate
    t += sim_dt;
    sim.update(planned_force);
  }

  // Final check
  double ref_pos = ref_pos_func(t);
  EXPECT_LT(std::abs(sim.state_[0] - ref_pos), 1e-2); // [m]
  EXPECT_LT(std::abs(sim.state_[1]), 1e-2); // [m/s]

  Eigen::Map<Eigen::VectorXd> computation_duration_vec(computation_duration_list.data(),
                                                       computation_duration_list.size());
  std::cout << "Computation time per control cycle:\n"
            << "  mean: " << computation_duration_vec.mean() << " [ms], stdev: "
            << std::sqrt((computation_duration_vec.array() - computation_duration_vec.mean()).square().mean())
            << " [ms], max: " << computation_duration_vec.maxCoeff() << " [ms]" << std::endl;

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:4 w lp # pos\n"
            << "  plot \"" << file_path << "\" u 1:5 w lp # force\n"
            << "  plot \"" << file_path << "\" u 1:6 w lp # computation_time\n";
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
