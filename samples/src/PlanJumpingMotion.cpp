/* Author: Masaki Murooka */

#include <fstream>

#include <CCC/Constants.h>
#include <CCC/LinearMpcTZ.h>

int main(int argc, char ** argv)
{
  double mass = 104.835; // [kg]
  double dt = 0.02; // [ms]
  double motion_duration = 5.0; // [s]

  int seq_len = static_cast<int>(motion_duration / dt);
  std::vector<bool> contact_seq(seq_len);
  Eigen::VectorXd ref_pos_seq(seq_len);
  for(int i = 0; i < seq_len; i++)
  {
    double t = i * dt;

    if(2.5 < t && t < 2.75)
    {
      contact_seq[i] = false;
    }
    else
    {
      contact_seq[i] = true;
    }

    ref_pos_seq[i] = 0.95;
  }
  double initial_pos = 0.95; // [m]
  double initial_vel = 0.0; // [m/s]

  CCC::LinearMpcTZ mpc(mass, dt);
  mpc.runLoop(contact_seq, initial_pos, initial_vel, ref_pos_seq);

  double start_time = 3.0; // [s]
  std::ofstream ofs_data("/tmp/PlanJumpingMotionData.txt");
  std::ofstream ofs_src("/tmp/PlanJumpingMotionSrc.txt");
  ofs_data << "time contact ref_pos planned_pos planned_vel planned_force" << std::endl;
  for(int i = 0; i < seq_len; i++)
  {
    double t = start_time + i * dt;
    ofs_data << t << " " << contact_seq[i] << " " << ref_pos_seq[i] << " " << mpc.planned_pos_seq_[i] << " "
             << mpc.planned_vel_seq_[i] << " " << mpc.planned_force_seq_[i] << std::endl;
    ofs_src << "{" << t << ", {{0, 0, 0}, {0, 0, " << mpc.planned_force_seq_[i] - mass * CCC::constants::g << "}}},"
            << std::endl;
  }
  std::cout << "Run the following commands in gnuplot:\n"
               "  set key autotitle columnhead\n"
               "  set key noenhanced\n"
               "  plot \"/tmp/PlanJumpingMotionData.txt\" u 1:3 w lp, \"\" u 1:4 w lp\n";

  return 0;
}
