/* Author: Masaki Murooka */

#pragma once

#include <CCC/VariantSequentialExtension.h>
#include <optional>

namespace CCC
{
/** \brief Linear MPC based on one-dimensional discrete dynamics on step switching.

    See the following for a detailed formulation.
      - S Xin, et al. Online relative footstep optimization for legged robots dynamic walking using discrete-time model
   predictive control. IROS, 2019.
 */
class StepMpc1d
{
  friend class StepMpc;

public:
  /** \brief Type of state-space model with fixed state dimension. */
  using _StateSpaceModel = StateSpaceModel<2, Eigen::Dynamic, Eigen::Dynamic>;

  /** \brief Discrete dynamics on step switching. */
  class StepModel : public _StateSpaceModel
  {
  public:
    /** \brief Constructor.
        \param com_height height of robot CoM [m]
        \param step_duration step duration [sec]
    */
    StepModel(double com_height, double step_duration);
  };

  /** \brief Reference data. */
  struct RefData
  {
    /** \brief Element of reference data. */
    struct Element
    {
      //! Whether it is in single support phase
      bool is_single_support = true;

      //! ZMP [m]
      double zmp = 0;

      //! End time [sec]
      double end_time = 0;
    };

    /** \brief List of reference element.

        The number of elements must be at least one.
        Consecutive double support phases are not allowed. (In contrast, consecutive single support phases are allowed.)
     */
    std::vector<Element> element_list;
  };

  /** \brief Planned data. */
  struct PlannedData
  {
    //! Current ZMP [m]
    double current_zmp = 0;

    /** \brief ZMP of next foot [m]

        null if the single support phase of the next foot is not included in the horizon
    */
    std::optional<double> next_foot_zmp;
  };

  /** \brief Initial parameter.

      First element is CoM position, and second element is CoM velocity.
  */
  using InitialParam = Eigen::Vector2d;

  /** \brief Weight parameter. */
  struct WeightParam
  {
    //! Weight of free ZMP (for future contacts)
    double free_zmp;

    //! Weight of fixed ZMP (for existing contacts)
    double fixed_zmp;

    //! Weight of ZMP in double support phase
    double double_support;

    //! Weight of CoM position
    double pos;

    //! Weight of CoM velocity
    double vel;

    //! Weight of absolute position of capture point
    double capture_point_abs;

    //! Weight of relative position of capture point and ZMP
    double capture_point_rel;

    /** \brief Constructor.
        \param _free_zmp weight of free ZMP (for future contacts)
        \param _fixed_zmp weight of fixed ZMP (for existing contacts)
        \param _double_support weight of ZMP in double support phase
        \param _pos weight of CoM position
        \param _vel weight of CoM velocity
        \param _capture_point_abs weight of absolute position of capture point
        \param _capture_point_rel weight of relative position of capture point and ZMP
    */
    WeightParam(double _free_zmp = 1e-2,
                double _fixed_zmp = 1e0,
                double _double_support = 1e0,
                double _pos = 0.0,
                double _vel = 0.0,
                double _capture_point_abs = 1e1,
                double _capture_point_rel = 1e1)
    : free_zmp(_free_zmp), fixed_zmp(_fixed_zmp), double_support(_double_support), pos(_pos), vel(_vel),
      capture_point_abs(_capture_point_abs), capture_point_rel(_capture_point_rel)
    {
    }
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param weight_param objective weight parameter
   */
  StepMpc1d(double com_height, const WeightParam & weight_param = WeightParam())
  : com_height_(com_height), weight_param_(weight_param)
  {
  }

  /** \brief Plan one step.
      \param ref_data reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned ZMP
   */
  PlannedData planOnce(const RefData & ref_data, const InitialParam & initial_param, double current_time);

protected:
  //! Height of robot CoM [m]
  double com_height_ = 0;

  //! Weight parameter
  WeightParam weight_param_;

  //! Sequential extension of state-space model
  std::shared_ptr<VariantSequentialExtension<2>> seq_ext_;
};

/** \brief Linear MPC based on discrete dynamics on step switching.

    See the following for a detailed formulation.
      - S Xin, et al. Online relative footstep optimization for legged robots dynamic walking using discrete-time model
   predictive control. IROS, 2019.

   \todo It is assumed that the left and right feet alternate in stepping (i.e., the same foot does not step in
   succession).
 */
class StepMpc
{
public:
  /** \brief Reference data. */
  struct RefData
  {
    /** \brief Element of reference data. */
    struct Element
    {
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      //! Whether it is in single support phase
      bool is_single_support = true;

      //! ZMP [m]
      Eigen::Vector2d zmp = Eigen::Vector2d::Zero();

      //! End time [sec]
      double end_time = 0;
    };

    /** \brief List of reference element.

        The number of elements must be at least one.
        Consecutive double support phases are not allowed. (In contrast, consecutive single support phases are allowed.)
     */
    std::vector<Element> element_list;
  };

  /** \brief Planned data. */
  struct PlannedData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Current ZMP [m]
    Eigen::Vector2d current_zmp = Eigen::Vector2d::Zero();

    /** \brief ZMP of next foot [m]

        null if the single support phase of the next foot is not included in the horizon
    */
    std::optional<Eigen::Vector2d> next_foot_zmp;
  };

  /** \brief Initial parameter. */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! CoM position [m]
    Eigen::Vector2d pos = Eigen::Vector2d::Zero();

    //! CoM velocity [m/s]
    Eigen::Vector2d vel = Eigen::Vector2d::Zero();
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param weight_param objective weight parameter
   */
  StepMpc(double com_height, const StepMpc1d::WeightParam & weight_param = StepMpc1d::WeightParam())
  : mpc_1d_(std::make_shared<StepMpc1d>(com_height, weight_param))
  {
  }

  /** \brief Plan one step.
      \param ref_data reference data
      \param initial_param initial parameter
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned ZMP
   */
  PlannedData planOnce(const RefData & ref_data, const InitialParam & initial_param, double current_time);

protected:
  //! One-dimensional linear MPC
  std::shared_ptr<StepMpc1d> mpc_1d_;
};
} // namespace CCC
