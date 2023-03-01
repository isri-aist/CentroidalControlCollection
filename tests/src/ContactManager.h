/* Author: Masaki Murooka */

#pragma once

#include <ForceColl/Contact.h>

/** \brief Make contact list from rectangle support region.
    \param rect_min_max min/max vertices of rectangle support region
*/
inline std::shared_ptr<ForceColl::Contact> makeContactFromRect(const std::array<Eigen::Vector2d, 2> & rect_min_max)
{
  std::vector<Eigen::Vector3d> vertex_list(4);
  vertex_list[0] << rect_min_max[0], 0.0;
  vertex_list[1] << rect_min_max[0][0], rect_min_max[1][1], 0.0;
  vertex_list[2] << rect_min_max[1], 0.0;
  vertex_list[3] << rect_min_max[1][0], rect_min_max[0][1], 0.0;

  constexpr double fricCoeff = 0.5;
  return std::make_shared<ForceColl::SurfaceContact>("ContactFromRect", fricCoeff, vertex_list,
                                                     sva::PTransformd::Identity());
}
