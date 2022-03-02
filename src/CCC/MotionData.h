/* Author: Masaki Murooka */

#pragma once

namespace CCC
{
/** \brief Motion data. */
template<class PosType, class VelType, class ForceType>
struct MotionDataBase
{
  //! Time [s]
  double time;

  //! Reference position
  PosType ref_pos;

  //! Reference velocity
  VelType ref_vel;

  //! Reference acceleration
  VelType ref_acc;

  //! Planned position
  PosType planned_pos;

  //! Planned velocity
  VelType planned_vel;

  //! Planned acceleration
  VelType planned_acc;

  //! Planned force
  ForceType planned_force;
};
} // namespace CCC
