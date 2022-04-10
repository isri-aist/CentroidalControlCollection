/* Author: Masaki Murooka */

#pragma once

#include <CCC/StateSpaceModel.h>

namespace CCC
{
/** \brief State-space model of CoM-ZMP dynamics with CoM jerk input and ZMP output.

      Dynamics is expressed by the following equation.
      \f{align*}{
      \ddot{c}_x = \dfrac{g}{c_z} (c_x - z_x)
      \f}
      \f$\boldsymbol{c}\f$ and \f$\boldsymbol{z}\f$ are CoM and ZMP.

      This can be represented as a linear time-invariant system as follows.
      \f{align*}{
      \boldsymbol{\dot{x}} &=
      \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix}
      \boldsymbol{x} +
      \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} u \\
      y &= \begin{bmatrix} 1 & 0 & - \dfrac{c_z}{g} \end{bmatrix} \boldsymbol{x}
      \f}

      State, control input, and output are expressed as follows.
      \f{align*}{
      \boldsymbol{x} = \begin{bmatrix} c_x \\ \dot{c}_x \\ \ddot{c}_x \end{bmatrix},
      u = \dddot{c}_x, y = z_x
      \f}
   */
class ComZmpModelJerkInput : public StateSpaceModel<3, 1, 1>
{
public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
  */
  ComZmpModelJerkInput(double com_height);
};
} // namespace CCC
