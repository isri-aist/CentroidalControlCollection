/* Author: Masaki Murooka */

#include <CCC/DcmTracking.h>

using namespace CCC;

Eigen::Vector2d DcmTracking::planOnce(const RefData & ref_data,
                                      const InitialParam & initial_param,
                                      double current_time) const
{
  // Calculate DCM at switching
  Eigen::Vector2d dcm_switch;
  if(ref_data.time_zmp_list.empty())
  {
    dcm_switch = ref_data.current_zmp;
  }
  else
  {
    // Check time
    for(const auto & time_zmp_kv : ref_data.time_zmp_list)
    {
      if(time_zmp_kv.first < current_time)
      {
        throw std::runtime_error("ZMP switching time must be in the future: " + std::to_string(time_zmp_kv.first)
                                 + " < " + std::to_string(current_time));
      }
    }

    dcm_switch = ref_data.time_zmp_list.rbegin()->second;
    for(auto time_zmp_it = std::next(ref_data.time_zmp_list.rbegin()); time_zmp_it != ref_data.time_zmp_list.rend();
        time_zmp_it++)
    {
      double zmp_duration = std::prev(time_zmp_it)->first - time_zmp_it->first;
      // Equation (18) in the paper
      dcm_switch = time_zmp_it->second + std::exp(-1 * omega_ * zmp_duration) * (dcm_switch - time_zmp_it->second);
    }
  }

  // Calculate target DCM: equation (19) in the paper
  Eigen::Vector2d target_dcm =
      ref_data.current_zmp
      + std::exp(omega_ * (current_time - ref_data.time_zmp_list.begin()->first)) * (dcm_switch - ref_data.current_zmp);

  // Calculate control ZMP: equation (24) in the paper
  Eigen::Vector2d control_zmp = ref_data.current_zmp + (1.0 + feedback_gain_ / omega_) * (initial_param - target_dcm);

  return control_zmp;
}
