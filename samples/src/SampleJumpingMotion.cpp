/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/LinearMpcXY.h>
#include <CCC/LinearMpcZ.h>

std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> makeVertexRidgeListFromRect(
    const std::array<Eigen::Vector2d, 2> & rect_min_max)
{
  std::vector<Eigen::Vector3d> vertexList(4);
  vertexList[0] << rect_min_max[0], 0.0;
  vertexList[1] << rect_min_max[0][0], rect_min_max[1][1], 0.0;
  vertexList[2] << rect_min_max[1], 0.0;
  vertexList[3] << rect_min_max[1][0], rect_min_max[0][1], 0.0;

  std::vector<Eigen::Vector3d> ridgeList(4);
  for(int i = 0; i < 4; i++)
  {
    double theta = 2 * M_PI * (static_cast<double>(i) / 4);
    ridgeList[i] << 0.5 * std::cos(theta), 0.5 * std::sin(theta), 1;
    ridgeList[i].normalize();
  }

  std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> vertex_ridge_list;
  for(const auto & vertex : vertexList)
  {
    for(const auto & ridge : ridgeList)
    {
      vertex_ridge_list.emplace_back(vertex, ridge);
    }
  }

  return vertex_ridge_list;
}

int main(int argc, char ** argv)
{
  double mass = 100.0; // [kg]
  double sim_dt = 0.03; // [sec]
  std::pair<double, double> motion_time_range(0.0, 3.0); // ([sec], [sec])

  double jump_start_t = 1.4; // [sec]
  double jump_end_t = 1.6; // [sec]
  double ref_com_z = 1.0; // [m]

  // Plan Z
  std::shared_ptr<CCC::LinearMpcZ> mpc_z;
  {
    std::function<bool(double)> contact_func = [jump_start_t, jump_end_t](double t) {
      return !(jump_start_t < t && t < jump_end_t);
    };
    std::function<double(double)> ref_pos_func = [ref_com_z](double t) { return ref_com_z; }; // [m]
    CCC::LinearMpcZ::InitialParam initial_pos_vel(ref_com_z, 0.0); // ([m], [m/s])
    double horizon_dt = 0.03; // [sec]
    double horizon_duration = 3.0; // [sec]

    mpc_z = std::make_shared<CCC::LinearMpcZ>(mass, horizon_dt);
    mpc_z->planLoop(contact_func, ref_pos_func, initial_pos_vel, motion_time_range, horizon_duration, sim_dt);

    mpc_z->dumpMotionDataSeq("/tmp/SampleJumpingMotionZ.txt", true);
  }

  // Plan XY
  std::shared_ptr<CCC::LinearMpcXY> mpc_xy;
  {
    std::function<CCC::LinearMpcXY::MotionParam(double)> motion_param_func =
        [mass, ref_com_z, jump_start_t, jump_end_t, sim_dt, motion_time_range, mpc_z](double t) {
          CCC::LinearMpcXY::MotionParam motion_param;

          // Set com_z and total_force_z
          int motion_data_z_idx = std::min(static_cast<int>(std::round((t - motion_time_range.first) / sim_dt)),
                                           static_cast<int>(mpc_z->motion_data_seq_.size()) - 1);
          const auto & motion_data_z = mpc_z->motion_data_seq_[motion_data_z_idx];
          motion_param.com_z = motion_data_z.planned_pos;
          motion_param.total_force_z = motion_data_z.planned_force;
          const auto & last_motion_data_z = mpc_z->motion_data_seq_[mpc_z->motion_data_seq_.size() - 1];
          if(std::abs(motion_data_z.time - std::min(t, last_motion_data_z.time)) > 1e-10)
          {
            throw std::runtime_error("motion_data_z.time is inconsistent. " + std::to_string(motion_data_z.time)
                                     + " != " + std::to_string(std::min(t, last_motion_data_z.time)));
          }

          // Set vertex_ridge_list
          if(t <= jump_start_t)
          {
            motion_param.vertex_ridge_list =
                makeVertexRidgeListFromRect({Eigen::Vector2d(-0.1, -0.1), Eigen::Vector2d(0.1, 0.1)}); // [m]
          }
          else if(t >= jump_end_t)
          {
            motion_param.vertex_ridge_list =
                makeVertexRidgeListFromRect({Eigen::Vector2d(0.4, -0.1), Eigen::Vector2d(0.6, 0.1)});
          }

          return motion_param;
        };
    std::function<CCC::LinearMpcXY::RefData(double)> ref_data_func = [jump_start_t, jump_end_t](double t) {
      CCC::LinearMpcXY::RefData ref_data;
      Eigen::Vector2d start_pos(0.0, 0.0); // [m]
      Eigen::Vector2d end_pos(0.5, 0.0); // [m]
      ref_data.pos =
          (end_pos - start_pos) / (jump_end_t - jump_start_t) * std::min(std::max(t, jump_start_t), jump_end_t)
          + (jump_end_t * start_pos - jump_start_t * end_pos) / (jump_end_t - jump_start_t);
      return ref_data;
    };
    CCC::LinearMpcXY::InitialParam initial_param;
    initial_param.pos << 0.0, 0.0; // [m]
    double horizon_dt = 0.03; // [sec]
    double horizon_duration = 1.0; // [sec]

    mpc_xy = std::make_shared<CCC::LinearMpcXY>(mass, horizon_dt);
    mpc_xy->planLoop(motion_param_func, ref_data_func, initial_param, motion_time_range, horizon_duration, sim_dt);

    mpc_xy->dumpMotionDataSeq("/tmp/SampleJumpingMotionXY.txt", true);
  }

  return 0;
}
