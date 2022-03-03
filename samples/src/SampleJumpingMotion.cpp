/* Author: Masaki Murooka */

#include <CCC/LinearMpcZ.h>

int main(int argc, char ** argv)
{
  double mass = 104.835; // [kg]
  double horizon_dt = 0.02; // [ms]

  CCC::LinearMpcZ mpc(mass, horizon_dt);

  std::function<bool(double)> contact_func = [](double t) { return !(5.5 < t && t < 5.7); };
  std::function<double(double)> ref_pos_func = [](double t) { return 0.85; }; // [m]
  Eigen::Vector2d initial_pos_vel(0.95, 0.0); // ([m], [m/s])
  std::pair<double, double> motion_time_range(3.0, 8.0); // ([s], [s])
  double horizon_duration = 4.0; // [sec]
  double sim_dt = 0.02; // [ms]

  mpc.planLoop(contact_func, ref_pos_func, initial_pos_vel, motion_time_range, horizon_duration, sim_dt);

  mpc.dumpMotionDataSeq("/tmp/PlanJumpingMotion.txt", true);

  return 0;
}
