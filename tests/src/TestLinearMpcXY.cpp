/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <CCC/Constants.h>
#include <CCC/LinearMpcXY.h>

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

TEST(TestLinearMpcXY, Test1)
{
  double mass = 100.0; // [kg]
  double horizon_dt = 0.1; // [sec]

  CCC::LinearMpcXY mpc(mass, horizon_dt);

  std::function<CCC::LinearMpcXY::MotionParam(double)> motion_param_func = [mass](double t) {
    CCC::LinearMpcXY::MotionParam motion_param;
    motion_param.com_z = 1.0; // [m]
    motion_param.total_force_z = mass * CCC::constants::g; // [N]
    std::array<Eigen::Vector2d, 2> rect_min_max;
    if(t < 3.0)
    {
      rect_min_max = {Eigen::Vector2d(0.9, -0.15), Eigen::Vector2d(1.1, 0.15)}; // [m]
    }
    else if(t < 4.0)
    {
      rect_min_max = {Eigen::Vector2d(0.9, 0.05), Eigen::Vector2d(1.1, 0.15)};
    }
    else if(t < 5.0)
    {
      rect_min_max = {Eigen::Vector2d(1.15, -0.15), Eigen::Vector2d(1.35, -0.05)};
    }
    else if(t < 6.0)
    {
      rect_min_max = {Eigen::Vector2d(1.4, 0.05), Eigen::Vector2d(1.6, 0.15)};
    }
    else
    {
      rect_min_max = {Eigen::Vector2d(1.4, -0.15), Eigen::Vector2d(1.6, 0.15)};
    }
    motion_param.vertex_ridge_list = makeVertexRidgeListFromRect(rect_min_max);
    return motion_param;
  };
  std::function<CCC::LinearMpcXY::RefData(double)> ref_data_func = [](double t) {
    CCC::LinearMpcXY::RefData ref_data;
    if(t < 3.0)
    {
      ref_data.pos << 1.0, 0.0; // [m]
    }
    else if(t < 4.0)
    {
      ref_data.pos << 1.0, 0.1;
    }
    else if(t < 5.0)
    {
      ref_data.pos << 1.25, -0.1;
    }
    else if(t < 6.0)
    {
      ref_data.pos << 1.5, 0.1;
    }
    else
    {
      ref_data.pos << 1.5, 0.0;
    }
    return ref_data;
  };
  CCC::LinearMpcXY::InitialParam initial_param;
  initial_param.pos << 1.0, 0.0; // [m]
  std::pair<double, double> motion_time_range(0.0, 8.0); // ([sec], [sec])
  double horizon_duration = 1.5; // [sec]
  double sim_dt = 0.05; // [sec]

  mpc.planLoop(motion_param_func, ref_data_func, initial_param, motion_time_range, horizon_duration, sim_dt);

  mpc.dumpMotionDataSeq("/tmp/TestLinearMpcXY.txt", true);

  // Check final state
  const auto & motion_data_final = mpc.motionDataSeq().rbegin()->second;
  EXPECT_LT((motion_data_final.planned_pos - motion_data_final.ref_pos).norm(),
            1e-3); // [m]
  EXPECT_LT((motion_data_final.planned_vel - motion_data_final.ref_vel).norm(),
            1e-3); // [m/s]
  EXPECT_LT((motion_data_final.planned_angular_momentum - motion_data_final.ref_angular_momentum).norm(),
            1e-3); // [kg m^2/s]
  const auto & ref_data_final = ref_data_func(motion_time_range.second);
  EXPECT_LT((motion_data_final.ref_pos - ref_data_final.pos).norm(),
            1e-10); // [m]
  EXPECT_LT((motion_data_final.ref_vel - ref_data_final.vel).norm(),
            1e-10); // [m/s]
  EXPECT_LT((motion_data_final.ref_angular_momentum - ref_data_final.angular_momentum).norm(),
            1e-10); // [kg m^2/s]
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
