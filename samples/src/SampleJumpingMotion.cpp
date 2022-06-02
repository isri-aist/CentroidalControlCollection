/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/LinearMpcXY.h>
#include <CCC/LinearMpcZ.h>

Eigen::Matrix<double, 6, Eigen::Dynamic> makeVertexRidgeListFromRect(const std::array<Eigen::Vector2d, 2> & rect_min_max)
{
  std::vector<Eigen::Vector3d> vertex_list(4);
  vertex_list[0] << rect_min_max[0], 0.0;
  vertex_list[1] << rect_min_max[0][0], rect_min_max[1][1], 0.0;
  vertex_list[2] << rect_min_max[1], 0.0;
  vertex_list[3] << rect_min_max[1][0], rect_min_max[0][1], 0.0;

  std::vector<Eigen::Vector3d> ridge_list(4);
  for(int i = 0; i < 4; i++)
  {
    double theta = 2 * M_PI * (static_cast<double>(i) / 4);
    ridge_list[i] << 0.5 * std::cos(theta), 0.5 * std::sin(theta), 1;
    ridge_list[i].normalize();
  }

  Eigen::Matrix<double, 6, Eigen::Dynamic> vertex_ridge_list(6, vertex_list.size() * ridge_list.size());
  int col_idx = 0;
  for(const auto & vertex : vertex_list)
  {
    for(const auto & ridge : ridge_list)
    {
      vertex_ridge_list.col(col_idx) << vertex, ridge;
      col_idx++;
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
    std::function<CCC::LinearMpcXY::MotionParam(double)> motion_param_func = [jump_start_t, jump_end_t,
                                                                              mpc_z](double t) {
      CCC::LinearMpcXY::MotionParam motion_param;

      // Set com_z and total_force_z
      {
        auto it = mpc_z->motionDataSeq().upper_bound(t);
        if(it == mpc_z->motionDataSeq().begin())
        {
          motion_param.com_z = it->second.planned_pos;
          motion_param.total_force_z = it->second.planned_force;
        }
        else if(it == mpc_z->motionDataSeq().end())
        {
          auto end_it = mpc_z->motionDataSeq().rbegin();
          motion_param.com_z = end_it->second.planned_pos;
          motion_param.total_force_z = end_it->second.planned_force;
        }
        else
        {
          double end_time = it->first;
          const auto & end_motion_data = it->second;
          it--;
          double start_time = it->first;
          const auto & start_motion_data = it->second;
          double ratio = (t - start_time) / (end_time - start_time);
          motion_param.com_z = (1 - ratio) * start_motion_data.planned_pos + ratio * end_motion_data.planned_pos;
          motion_param.total_force_z =
              (1 - ratio) * start_motion_data.planned_force + ratio * end_motion_data.planned_force;
        }
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
    double horizon_dt = 0.04; // [sec]
    double horizon_duration = 1.0; // [sec]

    mpc_xy = std::make_shared<CCC::LinearMpcXY>(mass, horizon_dt);
    mpc_xy->planLoop(motion_param_func, ref_data_func, initial_param, motion_time_range, horizon_duration, sim_dt);

    mpc_xy->dumpMotionDataSeq("/tmp/SampleJumpingMotionXY.txt", true);
  }

  return 0;
}
