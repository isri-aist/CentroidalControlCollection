/* Author: Masaki Murooka */

#pragma once

#include <deque>
#include <map>
#include <memory>
#include <unordered_map>

#include <Eigen/Core>

#include <ros/console.h>

#include <CCC/DcmTracking.h>
#include <CCC/FootGuidedControl.h>
#include <CCC/IntrinsicallyStableMpc.h>
#include <CCC/LinearMpcZmp.h>

/** \brief Foot. */
enum class Foot
{
  //! Left foot
  Left = 0,

  //! Right foot
  Right
};

/** \brief Get opposite foot. */
inline Foot opposite(const Foot & foot)
{
  return foot == Foot::Left ? Foot::Right : Foot::Left;
}

/** \brief Footstep. */
struct Footstep
{
  //! Foot
  Foot foot = Foot::Left;

  //! Foot position
  Eigen::Vector2d pos = Eigen::Vector2d::Zero();

  //! Time to start ZMP transition to support foot
  double transit_start_time = 0;

  //! Time to start foot swing
  double swing_start_time = 0;

  //! Time to end foot swing
  double swing_end_time = 0;

  //! Time to end ZMP transition to stance center
  double transit_end_time = 0;

  /** \brief Constructor. */
  Footstep() {}

  /** \brief Constructor.
      \param _foot foot
      \param _pos foot position
      \param _transit_start_time time to start to move ZMP to support foot
      \param transit_duration ZMP transition duration
      \param swing_duration foot swing duration
   */
  Footstep(const Foot & _foot,
           const Eigen::Vector2d & _pos,
           double _transit_start_time,
           double transit_duration,
           double swing_duration)
  : transit_start_time(_transit_start_time), swing_start_time(_transit_start_time + 0.5 * transit_duration),
    swing_end_time(_transit_start_time + 0.5 * transit_duration + swing_duration),
    transit_end_time(_transit_start_time + transit_duration + swing_duration), foot(_foot), pos(_pos)
  {
  }
};

/** \brief Footstance.
 *
 *  Stance means a set of contacts that are active at the same time.
 */
struct Footstance : public std::unordered_map<Foot, Eigen::Vector2d>
{
  using std::unordered_map<Foot, Eigen::Vector2d>::unordered_map;

  /** \brief Get stance center. */
  inline Eigen::Vector2d midPos() const
  {
    if(this->size() == 0)
    {
      ROS_WARN("[Footstance::midPos] The number of contacts is zero.");
      return Eigen::Vector2d::Zero();
    }
    else if(this->size() == 1)
    {
      return this->begin()->second;
    }
    else
    {
      return 0.5 * (this->at(Foot::Left) + this->at(Foot::Right));
    }
  }

  /** \brief Get min/max vertices of support region. */
  inline std::array<Eigen::Vector2d, 2> supportRegion() const
  {
    std::array<Eigen::Vector2d, 2> region_min_max;
    if(this->size() == 0)
    {
      ROS_WARN("[Footstance::supportRegion] The number of contacts is zero.");
      region_min_max[0] = Eigen::Vector2d::Zero();
      region_min_max[1] = Eigen::Vector2d::Zero();
    }
    else if(this->size() == 1)
    {
      region_min_max[0] = this->begin()->second;
      region_min_max[1] = this->begin()->second;
    }
    else
    {
      region_min_max[0] = this->at(Foot::Left).cwiseMin(this->at(Foot::Right));
      region_min_max[1] = this->at(Foot::Left).cwiseMax(this->at(Foot::Right));
    }
    return region_min_max;
  }
};

/** \brief Footstep manager. */
class FootstepManager
{
public:
  /** \brief Constructor.
      \param initial_footstance initial footstance
  */
  FootstepManager(const Footstance & initial_footstance = {{Foot::Left, Eigen::Vector2d(0.0, 0.1)},
                                                           {Foot::Right, Eigen::Vector2d(0.0, -0.1)}})
  : footstance_(initial_footstance)
  {
  }

  /** \brief Update.
      \param current_time current time

      This method should be called once every control cycle.
  */
  inline void update(double current_time)
  {
    // Update footstance
    if(!footstep_list_.empty() && footstep_list_.front().swing_end_time <= current_time)
    {
      footstance_.at(footstep_list_.front().foot) = footstep_list_.front().pos;
    }

    // Remove old footsteps
    while(!footstep_list_.empty() && footstep_list_.front().transit_end_time < current_time)
    {
      prev_footstep_ = std::make_shared<Footstep>(footstep_list_.front());
      footstep_list_.pop_front();
    }

    // Set ref_zmp_list and ref_footstance_list_
    ref_zmp_list_.clear();
    ref_footstance_list_.clear();
    if(footstep_list_.empty())
    {
      // If there is no footstep, keep current footstance
      ref_zmp_list_.emplace(current_time, footstance_.midPos());
      ref_zmp_list_.emplace(current_time + horizon_duration_, footstance_.midPos());
      ref_footstance_list_.emplace(current_time, footstance_);
      ref_footstance_list_.emplace(current_time + horizon_duration_, footstance_);
    }
    else
    {
      // If first footstep is not executed yet, set the current footstance first
      if(current_time < footstep_list_.front().transit_start_time)
      {
        ref_zmp_list_.emplace(current_time, footstance_.midPos());
        ref_footstance_list_.emplace(current_time, footstance_);
      }

      // Apply footsteps in horizon in order
      Footstance tmp_footstance = footstance_;
      for(auto footstep_it = footstep_list_.begin();
          footstep_it != footstep_list_.end() && footstep_it->transit_start_time <= current_time + horizon_duration_;
          footstep_it++)
      {
        ref_zmp_list_.emplace(footstep_it->transit_start_time, tmp_footstance.midPos());
        ref_footstance_list_.emplace(footstep_it->transit_start_time, tmp_footstance);

        tmp_footstance.erase(footstep_it->foot);
        ref_zmp_list_.emplace(footstep_it->swing_start_time, tmp_footstance.at(opposite(footstep_it->foot)));
        ref_footstance_list_.emplace(footstep_it->swing_start_time, tmp_footstance);

        tmp_footstance.emplace(footstep_it->foot, footstep_it->pos);
        ref_zmp_list_.emplace(footstep_it->swing_end_time, tmp_footstance.at(opposite(footstep_it->foot)));
        ref_footstance_list_.emplace(footstep_it->swing_end_time, tmp_footstance);

        ref_zmp_list_.emplace(footstep_it->transit_end_time, tmp_footstance.midPos());
      }

      // If last footstep is before the end of horizon, set the last footstance
      if(ref_zmp_list_.rbegin()->first < current_time + horizon_duration_)
      {
        ref_zmp_list_.emplace(current_time + horizon_duration_, tmp_footstance.midPos());
        ref_footstance_list_.emplace(current_time + horizon_duration_, tmp_footstance);
      }
    }
  }

  /** \brief Append footstep.
      \param footstep footstep
  */
  inline void appendFootstep(const Footstep & footstep)
  {
    if(!footstep_list_.empty() && footstep.transit_start_time < footstep_list_.back().transit_end_time)
    {
      throw std::runtime_error("[FootstepManager::appendFootstep] transit_start_time of specified footstep must be "
                               "after transit_end_time of last footstep.");
    }

    footstep_list_.push_back(footstep);
  }

  /** \brief Get reference ZMP.
      \param t time
  */
  inline Eigen::Vector2d refZmp(double t) const
  {
    // Add small values to avoid numerical instability at inequality bounds
    constexpr double epsilon_t = 1e-6;
    t += epsilon_t;

    auto ref_zmp_it = ref_zmp_list_.upper_bound(t);
    double ratio = (t - std::prev(ref_zmp_it)->first) / (ref_zmp_it->first - std::prev(ref_zmp_it)->first);
    return (1 - ratio) * std::prev(ref_zmp_it)->second + ratio * ref_zmp_it->second;
  }

  /** \brief Get ZMP limits.
      \param t time
  */
  inline std::array<Eigen::Vector2d, 2> zmpLimits(double t) const
  {
    const auto & footstance = std::prev(ref_footstance_list_.upper_bound(t))->second;
    std::array<Eigen::Vector2d, 2> zmp_limits;
    zmp_limits = footstance.supportRegion();
    zmp_limits[0] -= 0.5 * foot_size_;
    zmp_limits[1] += 0.5 * foot_size_;
    return zmp_limits;
  }

  /** \brief Make DcmTracking::RefData instance.
      \param current_time current time
      \param horizon_duration horizon duration
  */
  inline CCC::DcmTracking::RefData makeDcmTrackingRefData(double current_time, int horizon_duration = 5.0) const
  {
    CCC::DcmTracking::RefData ref_data;
    ref_data.current_zmp = std::prev(ref_zmp_list_.upper_bound(current_time))->second;
    for(auto time_zmp_it = ref_zmp_list_.upper_bound(current_time);
        time_zmp_it != ref_zmp_list_.end() && time_zmp_it->first < current_time + horizon_duration; time_zmp_it++)
    {
      ref_data.time_zmp_list.emplace(*time_zmp_it);
    }
    return ref_data;
  }

  /** \brief Make FootGuidedControl::RefData instance.
      \param current_time current time
  */
  inline CCC::FootGuidedControl::RefData makeFootGuidedControlRefData(double current_time) const
  {
    double constant_zmp_duration = 1.0; // [sec]
    double footsteps_concat_duration_thre = 0.1; // [sec]
    double horizon_margin = 1e-3; // [sec]
    CCC::FootGuidedControl::RefData ref_data;

    if(footstep_list_.empty())
    {
      ref_data.transit_start_zmp = footstance_.midPos();
      ref_data.transit_end_zmp = ref_data.transit_start_zmp;
      ref_data.transit_start_time = current_time + constant_zmp_duration;
      ref_data.transit_duration = 0;
    }
    else
    {
      const auto & footstep = footstep_list_.front();
      if(current_time < footstep.swing_start_time)
      {
        ref_data.transit_start_zmp = footstance_.midPos();
        ref_data.transit_end_zmp = footstance_.at(opposite(footstep.foot));
        ref_data.transit_start_time = footstep.transit_start_time;
        ref_data.transit_duration = footstep.swing_start_time - footstep.transit_start_time;

        // If the double support duration is short, concatenate the previous and current footsteps to avoid the horizon
        // becoming too short
        if(prev_footstep_)
        {
          if(prev_footstep_->foot == opposite(footstep.foot)
             && footstep.transit_start_time - prev_footstep_->transit_end_time < footsteps_concat_duration_thre)
          {
            ref_data.transit_start_zmp = footstance_.at(opposite(prev_footstep_->foot));
            ref_data.transit_start_time = prev_footstep_->swing_end_time;
            ref_data.transit_duration = footstep.swing_start_time - prev_footstep_->swing_end_time;
          }
        }
      }
      else
      {
        ref_data.transit_start_zmp = footstance_.at(opposite(footstep.foot));
        Footstance tmp_footstance = footstance_;
        tmp_footstance.at(footstep.foot) = footstep.pos;
        ref_data.transit_end_zmp = tmp_footstance.midPos();
        ref_data.transit_start_time = footstep.swing_end_time;
        ref_data.transit_duration = footstep.transit_end_time - footstep.swing_end_time;

        // If the double support duration is short, concatenate the current and next footsteps to avoid the horizon
        // becoming too short
        if(footstep_list_.size() >= 2)
        {
          Footstep next_footstep = footstep_list_[1];
          if(next_footstep.foot == opposite(footstep.foot)
             && next_footstep.transit_start_time - footstep.transit_end_time < footsteps_concat_duration_thre)
          {
            ref_data.transit_end_zmp = footstep.pos;
            ref_data.transit_duration = next_footstep.swing_start_time - footstep.swing_end_time;
          }
        }
      }

      // Ensure a horizon, since a horizon close to zero produces a very large input
      if(ref_data.transit_start_time + ref_data.transit_duration < current_time + horizon_margin)
      {
        ref_data.transit_duration += horizon_margin;
      }
    }

    return ref_data;
  }

  /** \brief Make LinearMpcZmp::RefData instance.
      \param t time
  */
  inline CCC::LinearMpcZmp::RefData makeLinearMpcZmpRefData(double t) const
  {
    CCC::LinearMpcZmp::RefData ref_data;
    ref_data.zmp_limits = zmpLimits(t);
    return ref_data;
  }

  /** \brief Make IntrinsicallyStableMpc::RefData instance.
      \param t time
  */
  inline CCC::IntrinsicallyStableMpc::RefData makeIntrinsicallyStableMpcRefData(double t) const
  {
    CCC::IntrinsicallyStableMpc::RefData ref_data;
    ref_data.zmp = refZmp(t);
    ref_data.zmp_limits = zmpLimits(t);
    return ref_data;
  }

public:
  //! Footstance
  Footstance footstance_;

  //! Footstep list
  std::deque<Footstep> footstep_list_;

  //! Previous footstep
  std::shared_ptr<Footstep> prev_footstep_;

  //! Horizon duration [sec]
  double horizon_duration_ = 10.0;

  //! Size of foot region [m]
  Eigen::Vector2d foot_size_ = Eigen::Vector2d(0.1, 0.05);

protected:
  //! Reference ZMP list
  std::map<double, Eigen::Vector2d> ref_zmp_list_;

  //! Reference footstance list
  std::map<double, Footstance> ref_footstance_list_;
};
