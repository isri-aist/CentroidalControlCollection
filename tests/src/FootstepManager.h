/* Author: Masaki Murooka */

#pragma once

#include <deque>
#include <map>
#include <unordered_map>

#include <Eigen/Core>

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
    return 0.5 * (this->at(Foot::Left) + this->at(Foot::Right));
  }
};

/** \brief Footstep manager. */
class FootstepManager
{
public:
  /** \brief Constructor.
      \param initial_footstance initial footstance
  */
  FootstepManager(const Footstance & initial_footstance = {{Foot::Left, Eigen::Vector2d(0.1, 0.1)},
                                                           {Foot::Right, Eigen::Vector2d(-0.1, -0.1)}})
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
      footstep_list_.pop_front();
    }

    // Set ref_zmp_list
    ref_zmp_list_.clear();
    if(footstep_list_.empty())
    {
      // If there is no footstep, keep current footstance
      ref_zmp_list_.emplace(current_time, footstance_.midPos());
      ref_zmp_list_.emplace(current_time + horizon_duration_, footstance_.midPos());
    }
    else
    {
      // If first footstep is not executed yet, set the current footstance first
      if(current_time < footstep_list_.front().transit_start_time)
      {
        ref_zmp_list_.emplace(current_time, footstance_.midPos());
      }

      // Apply footsteps in horizon in order
      Footstance tmp_footstance = footstance_;
      for(auto footstep_it = footstep_list_.begin();
          footstep_it != footstep_list_.end() && footstep_it->transit_start_time <= current_time + horizon_duration_;
          footstep_it++)
      {
        ref_zmp_list_.emplace(footstep_it->transit_start_time, tmp_footstance.midPos());
        ref_zmp_list_.emplace(footstep_it->swing_start_time, tmp_footstance.at(opposite(footstep_it->foot)));
        ref_zmp_list_.emplace(footstep_it->swing_end_time, tmp_footstance.at(opposite(footstep_it->foot)));
        tmp_footstance.at(footstep_it->foot) = footstep_it->pos;
        ref_zmp_list_.emplace(footstep_it->transit_end_time, tmp_footstance.midPos());
      }

      // If last footstep is before the end of horizon, set the last footstance
      if(ref_zmp_list_.rbegin()->first < current_time + horizon_duration_)
      {
        ref_zmp_list_.emplace(current_time + horizon_duration_, tmp_footstance.midPos());
      }
    }
  }

  /** \brief Get reference zmp.
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

  /** \brief Append footstep.
      \param footstep footstep
  */
  inline void appendFootstep(const Footstep & footstep)
  {
    if(!footstep_list_.empty() && footstep.transit_start_time < footstep_list_.back().transit_end_time)
    {
      throw std::runtime_error(
          "transit_start_time of specified footstep must be after transit_end_time of last footstep.");
    }

    footstep_list_.push_back(footstep);
  }

public:
  //! Footstance
  Footstance footstance_;

  //! Footstep list
  std::deque<Footstep> footstep_list_;

  //! Horizon duration
  double horizon_duration_ = 10.0; // [sec]

protected:
  //! Reference ZMP list
  std::map<double, Eigen::Vector2d> ref_zmp_list_;
};
