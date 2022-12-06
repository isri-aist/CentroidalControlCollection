/* Author: Masaki Murooka */

#pragma once

#include <Eigen/Core>

/** \brief Make list of vertex and ridge from rectangle support region.
    \param rect_min_max min/max vertices of rectangle support region
*/
inline Eigen::Matrix<double, 6, Eigen::Dynamic> makeVertexRidgeListFromRect(
    const std::array<Eigen::Vector2d, 2> & rect_min_max)
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
