/* Author: Masaki Murooka */

#include <CCC/DcmTracking.h>

using namespace CCC;

Eigen::Vector2d DcmTracking::planOnce(const std::map<double, Eigen::Vector2d> & time_zmp_list,
                                      const InitialParam & initial_param,
                                      double current_time) const
{
  // Calculate DCM at switching
  Eigen::Vector2d dcm_switch;
  if(time_zmp_list.empty())
  {
    dcm_switch = initial_param.ref_zmp;
  }
  else
  {
    // Check time
    for(const auto & time_zmp_kv : time_zmp_list)
    {
      if(time_zmp_kv.first < current_time)
      {
        throw std::runtime_error("ZMP switching time must be in the future: " + std::to_string(time_zmp_kv.first)
                                 + " < " + std::to_string(current_time));
      }
    }

    dcm_switch = time_zmp_list.rbegin()->second;
    for(auto time_zmp_it = std::next(time_zmp_list.rbegin()); time_zmp_it != time_zmp_list.rend(); time_zmp_it++)
    {
      double zmp_duration = std::prev(time_zmp_it)->first - time_zmp_it->first;
      // Equation (18) in the paper
      dcm_switch = time_zmp_it->second + std::exp(-1 * omega_ * zmp_duration) * (dcm_switch - time_zmp_it->second);
    }
  }

  // Calculate target DCM: equation (19) in the paper
  Eigen::Vector2d target_dcm =
      initial_param.ref_zmp
      + std::exp(omega_ * (current_time - time_zmp_list.begin()->first)) * (dcm_switch - initial_param.ref_zmp);

  // Calculate control ZMP: equation (24) in the paper
  Eigen::Vector2d control_zmp =
      initial_param.ref_zmp + (1.0 + feedback_gain_ / omega_) * (initial_param.current_dcm - target_dcm);

  return control_zmp;
}
