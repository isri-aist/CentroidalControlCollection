/* Author: Masaki Murooka */

#include <cmath>
#include <stdexcept>
#include <string>

#include <CCC/FootGuidedControl.h>

using namespace CCC;

double FootGuidedControl::planOnce(const RefData & ref_data,
                                   const InitialParam & initial_param,
                                   double current_time) const
{
  double capture_point = initial_param;

  if(!(ref_data.transit_duration >= 0))
  {
    throw std::runtime_error("Transition duration must be non-negative: " + std::to_string(ref_data.transit_duration));
  }

  double transit_end_time = ref_data.transit_start_time + ref_data.transit_duration;
  constexpr double future_margin_duration = 1e-6;
  if(!(transit_end_time >= current_time + future_margin_duration))
  {
    throw std::runtime_error("Transition end time must be in the future with some margin: "
                             + std::to_string(transit_end_time) + " < " + std::to_string(current_time) + " + "
                             + std::to_string(future_margin_duration));
  }

  double planned_zmp;
  if(ref_data.transit_duration == 0)
  {
    // Equation (7) of Kojio's paper
    planned_zmp = ref_data.transit_start_zmp
                  + 2
                        * ((capture_point - ref_data.transit_start_zmp)
                           - (ref_data.transit_end_zmp - ref_data.transit_start_zmp)
                                 * std::exp(-1 * omega_ * (ref_data.transit_start_time - current_time)))
                        / (1.0 - std::exp(-2 * omega_ * (ref_data.transit_start_time - current_time)));
  }
  else
  {
    double zmp_transit_vel = (ref_data.transit_end_zmp - ref_data.transit_start_zmp) / ref_data.transit_duration;
    if(current_time <= ref_data.transit_start_time)
    {
      // Equation (17) of Kojio's paper
      planned_zmp = ref_data.transit_start_zmp
                    + (2 * (capture_point - ref_data.transit_start_zmp)
                       + 2 * zmp_transit_vel / omega_
                             * (std::exp(-1 * omega_ * (transit_end_time - current_time))
                                - std::exp(-1 * omega_ * (ref_data.transit_start_time - current_time))))
                          / (1.0 - std::exp(-2 * omega_ * (transit_end_time - current_time)));
    }
    else
    {
      double current_ref_zmp =
          ref_data.transit_start_zmp + zmp_transit_vel * (current_time - ref_data.transit_start_time);
      // Equation (17) of Kojio's paper
      planned_zmp =
          current_ref_zmp
          + (2 * (capture_point - current_ref_zmp)
             + 2 * zmp_transit_vel / omega_ * (std::exp(-1 * omega_ * (transit_end_time - current_time)) - 1.0))
                / (1.0 - std::exp(-2 * omega_ * (transit_end_time - current_time)));
    }
  }

  return planned_zmp;
}
