/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/StepMpc.h>

using namespace CCC;

StepMpc1d::StepModel::StepModel(double com_height, double step_duration) : _StateSpaceModel(2, 1, 3)
{
  double omega = std::sqrt(constants::g / com_height);
  double exp_omega = std::exp(omega * step_duration);
  double exp_omega_inv = 1.0 / exp_omega;
  Ad_.resize(2, 2);
  Ad_ << 0.5 * (exp_omega + exp_omega_inv), 0.5 * (exp_omega - exp_omega_inv) / omega,
      0.5 * omega * (exp_omega - exp_omega_inv), 0.5 * (exp_omega + exp_omega_inv);
  Bd_.resize(2, 1);
  Bd_ << 1.0 - 0.5 * (exp_omega + exp_omega_inv), 0.5 * omega * (exp_omega_inv - exp_omega);
  Ed_.setZero(2);

  C_.resize(3, 2);
  C_ << 1.0, 0.0, 0.0, 1.0, 1.0, 1.0 / omega;
  D_.setZero(3, 1);
  F_.setZero(3);
}

StepMpc1d::PlannedData StepMpc1d::planOnce(const RefData & ref_data,
                                           const InitialParam & initial_param,
                                           double current_time)
{
  size_t num_var = ref_data.element_list.size();
  Eigen::MatrixXd eq_mat = Eigen::MatrixXd::Zero(num_var, num_var);
  Eigen::MatrixXd eq_vec = Eigen::VectorXd::Zero(num_var);

  // Setup sequential extension
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list;
  for(size_t i = 0; i < ref_data.element_list.size(); i++)
  {
    double step_duration;
    if(i == 0)
    {
      step_duration = ref_data.element_list[i].end_time - current_time;
    }
    else
    {
      step_duration = ref_data.element_list[i].end_time - ref_data.element_list[i - 1].end_time;
    }
    model_list.push_back(std::make_shared<StepModel>(com_height_, step_duration));
  }
  seq_ext_ = std::make_shared<VariantSequentialExtension<2>>(model_list, true);

  // Set coefficients for ZMP
  Eigen::VectorXd ref_zmp_vec(num_var);
  for(size_t i = 0; i < num_var; i++)
  {
    ref_zmp_vec[i] = ref_data.element_list[i].zmp;
  }
  Eigen::VectorXd eq_mat_diag_zmp = Eigen::VectorXd::Constant(num_var, weight_param_.free_zmp);
  eq_mat_diag_zmp[0] = weight_param_.fixed_zmp;
  if(num_var > 1 && !ref_data.element_list[0].is_single_support)
  {
    eq_mat_diag_zmp[1] = weight_param_.fixed_zmp;
  }
  eq_mat.diagonal() += eq_mat_diag_zmp;
  eq_vec += -1 * eq_mat_diag_zmp.cwiseProduct(ref_zmp_vec);

  // Set coefficients for double support
  Eigen::Matrix3d eq_mat_double_support;
  eq_mat_double_support << 1.0, -2.0, 1.0, -2.0, 4.0, -2.0, 1.0, -2.0, 1.0;
  eq_mat_double_support *= weight_param_.double_support;
  for(size_t i = 0; i < num_var; i++)
  {
    if(i != 0 && i != num_var - 1)
    {
      if(ref_data.element_list[i - 1].is_single_support && !ref_data.element_list[i].is_single_support
         && ref_data.element_list[i + 1].is_single_support)
      {
        eq_mat.block<3, 3>(i - 1, i - 1) += eq_mat_double_support;
      }
    }
  }

  // Set coefficients for CoM position
  if(weight_param_.pos > 0.0)
  {
    Eigen::VectorXd ref_pos_vec(num_var);
    Eigen::MatrixXd select_mat_pos = Eigen::MatrixXd::Zero(num_var, 3 * num_var);
    for(size_t i = 0; i < num_var; i++)
    {
      ref_pos_vec[i] = ref_data.element_list[std::min(i + 1, num_var - 1)].zmp;
      select_mat_pos(i, 3 * i) = 1.0;
    }
    Eigen::MatrixXd eq_mat_pos_sub = weight_param_.pos * seq_ext_->B_seq_.transpose() * select_mat_pos.transpose();
    eq_mat += eq_mat_pos_sub * select_mat_pos * seq_ext_->B_seq_;
    eq_vec += eq_mat_pos_sub * (select_mat_pos * seq_ext_->A_seq_ * initial_param - ref_pos_vec);
  }

  // Set coefficients for CoM velocity
  if(weight_param_.vel > 0.0)
  {
    Eigen::MatrixXd select_mat_vel = Eigen::MatrixXd::Zero(num_var, 3 * num_var);
    for(size_t i = 0; i < num_var; i++)
    {
      select_mat_vel(i, 3 * i + 1) = 1.0;
    }
    Eigen::MatrixXd eq_mat_vel_sub =
        weight_param_.vel * seq_ext_->B_seq_.transpose() * select_mat_vel.transpose() * select_mat_vel;
    eq_mat += eq_mat_vel_sub * seq_ext_->B_seq_;
    eq_vec += eq_mat_vel_sub * seq_ext_->A_seq_ * initial_param;
  }

  // Set coefficients for absolute position of capture point
  if(weight_param_.capture_point_abs > 0.0)
  {
    size_t num_future_single_support = 0;
    for(size_t i = 0; i < num_var; i++)
    {
      if(i >= 1 && ref_data.element_list[i].is_single_support)
      {
        num_future_single_support++;
      }
    }
    Eigen::VectorXd ref_cp_vec(num_var);
    Eigen::MatrixXd select_mat_cp = Eigen::MatrixXd::Zero(num_var, 3 * num_var);
    if(num_future_single_support == 0)
    {
      for(size_t i = 0; i < num_var; i++)
      {
        if(i >= 1 || !ref_data.element_list[i].is_single_support)
        {
          ref_cp_vec[i] = ref_data.element_list[i].zmp;
          select_mat_cp(i, 3 * i + 2) = 1.0;
        }
      }
    }
    else
    {
      for(size_t i = 0; i < num_var; i++)
      {
        if(i >= 1 && ref_data.element_list[i].is_single_support)
        {
          ref_cp_vec[i] = ref_data.element_list[i].zmp;
          select_mat_cp(i, 3 * (i - 1) + 2) = 1.0;
        }
      }
    }
    Eigen::MatrixXd eq_mat_cp_sub =
        weight_param_.capture_point_abs * seq_ext_->B_seq_.transpose() * select_mat_cp.transpose();
    eq_mat += eq_mat_cp_sub * select_mat_cp * seq_ext_->B_seq_;
    eq_vec += eq_mat_cp_sub * (select_mat_cp * seq_ext_->A_seq_ * initial_param - ref_cp_vec);
  }

  // Set coefficients for relative position of capture point and ZMP
  if(weight_param_.capture_point_rel > 0.0)
  {
    size_t num_future_single_support = 0;
    for(size_t i = 0; i < num_var; i++)
    {
      if(i >= 1 && ref_data.element_list[i].is_single_support)
      {
        num_future_single_support++;
      }
    }
    if(num_future_single_support >= 1)
    {
      Eigen::MatrixXd select_mat_cp = Eigen::MatrixXd::Zero(num_var, 3 * num_var);
      Eigen::MatrixXd select_mat_zmp = Eigen::MatrixXd::Zero(num_var, num_var);
      for(size_t i = 0; i < num_var; i++)
      {
        if(i >= 1 && ref_data.element_list[i].is_single_support)
        {
          select_mat_cp(i, 3 * (i - 1) + 2) = 1.0;
          select_mat_zmp(i, i) = 1.0;
        }
      }
      Eigen::MatrixXd eq_mat_cp_sub =
          weight_param_.capture_point_rel * (select_mat_cp * seq_ext_->B_seq_ - select_mat_zmp).transpose();
      eq_mat += eq_mat_cp_sub * (select_mat_cp * seq_ext_->B_seq_ - select_mat_zmp);
      eq_vec += eq_mat_cp_sub * select_mat_cp * seq_ext_->A_seq_ * initial_param;
    }
  }

  // Solve equation
  Eigen::VectorXd eq_sol = eq_mat.colPivHouseholderQr().solve(-1 * eq_vec);

  // Set planned data
  PlannedData planned_data;
  planned_data.current_zmp = eq_sol[0];
  if(num_var > 1)
  {
    for(size_t i = 1; i < num_var; i++)
    {
      if(ref_data.element_list[i].is_single_support)
      {
        planned_data.next_foot_zmp = eq_sol[i];
        break;
      }
    }
  }
  return planned_data;
}

StepMpc::PlannedData StepMpc::planOnce(const RefData & ref_data,
                                       const InitialParam & initial_param,
                                       double current_time)
{
  PlannedData planned_data;
  StepMpc1d::RefData ref_data_1d;
  StepMpc1d::InitialParam initial_param_1d;
  StepMpc1d::PlannedData planned_data_1d;

  ref_data_1d.element_list.clear();
  for(const auto & element : ref_data.element_list)
  {
    StepMpc1d::RefData::Element element_1d;
    element_1d.is_single_support = element.is_single_support;
    element_1d.zmp = element.zmp.x();
    element_1d.end_time = element.end_time;
    ref_data_1d.element_list.push_back(element_1d);
  }
  initial_param_1d << initial_param.pos.x(), initial_param.vel.x();
  planned_data_1d = mpc_1d_->planOnce(ref_data_1d, initial_param_1d, current_time);
  planned_data.current_zmp.x() = planned_data_1d.current_zmp;
  if(planned_data_1d.next_foot_zmp)
  {
    planned_data.next_foot_zmp = Eigen::Vector2d::Zero();
    planned_data.next_foot_zmp->x() = planned_data_1d.next_foot_zmp.value();
  }

  ref_data_1d.element_list.clear();
  for(const auto & element : ref_data.element_list)
  {
    StepMpc1d::RefData::Element element_1d;
    element_1d.is_single_support = element.is_single_support;
    element_1d.zmp = element.zmp.y();
    element_1d.end_time = element.end_time;
    ref_data_1d.element_list.push_back(element_1d);
  }
  initial_param_1d << initial_param.pos.y(), initial_param.vel.y();
  planned_data_1d = mpc_1d_->planOnce(ref_data_1d, initial_param_1d, current_time);
  planned_data.current_zmp.y() = planned_data_1d.current_zmp;
  if(planned_data.next_foot_zmp)
  {
    planned_data.next_foot_zmp->y() = planned_data_1d.next_foot_zmp.value();
  }

  return planned_data;
}
