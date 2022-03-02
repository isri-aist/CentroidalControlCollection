/* Author: Masaki Murooka */

#include <fstream>

#include <CCC/Constants.h>
#include <CCC/EigenTypes.h>
#include <CCC/LinearMpcXY.h>

using namespace CCC;

Eigen::Vector6d LinearMpcXY::InitialParam::toState(double mass) const
{
  Eigen::Vector6d state;
  state << mass * pos.x(), mass * vel.x(), mass * pos.y(), mass * vel.y(), angular_momentum;
  return state;
}

Eigen::Vector6d LinearMpcXY::RefData::toOutput(double mass) const
{
  Eigen::Vector6d output;
  output << mass * pos.x(), mass * vel.x(), mass * pos.y(), mass * vel.y(), angular_momentum;
  return output;
}

Eigen::VectorXd LinearMpcXY::WeightParam::inputWeight(int total_input_dim) const
{
  return Eigen::VectorXd::Constant(total_input_dim, force);
}

Eigen::VectorXd LinearMpcXY::WeightParam::outputWeight(size_t seq_len) const
{
  Eigen::Vector6d output_weight_one;
  output_weight_one << linear_momentum_integral.x(), linear_momentum.x(), linear_momentum_integral.y(),
      linear_momentum.y(), angular_momentum.x(), angular_momentum.y();

  Eigen::VectorXd output_weight(RefData::outputDim() * seq_len);
  for(size_t i = 0; i < seq_len; i++)
  {
    output_weight.segment<RefData::outputDim()>(i * RefData::outputDim()) << output_weight_one;
  }
  return output_weight;
}

LinearMpcXY::Model::Model(double mass, const MotionParam & motion_param, int output_dim)
: StateSpaceModel(LinearMpcXY::state_dim_, motion_param.vertex_ridge_list.size(), output_dim),
  motion_param_(motion_param)
{
  A_(0, 1) = 1;
  A_(2, 3) = 1;
  A_(4, 2) = -1 * motion_param_.total_force_z / mass;
  A_(5, 0) = motion_param_.total_force_z / mass;

  for(size_t i = 0; i < motion_param_.vertex_ridge_list.size(); i++)
  {
    const auto & vertex = motion_param_.vertex_ridge_list[i].first;
    const auto & ridge = motion_param_.vertex_ridge_list[i].second;
    B_.col(i) << 0, ridge.x(), 0, ridge.y(),
        -1 * (vertex.z() - motion_param_.com_z) * ridge.y() + vertex.y() * ridge.z(),
        (vertex.z() - motion_param_.com_z) * ridge.x() + -1 * vertex.x() * ridge.z();
  }
}

LinearMpcXY::SimModel::SimModel(double mass, const MotionParam & motion_param) : Model(mass, motion_param, 8)
{
  C_(0, 0) = 1 / mass;
  C_(1, 1) = 1 / mass;
  C_(3, 2) = 1 / mass;
  C_(4, 3) = 1 / mass;
  C_(6, 4) = 1;
  C_(7, 5) = 1;

  for(size_t i = 0; i < motion_param_.vertex_ridge_list.size(); i++)
  {
    const auto & ridge = motion_param_.vertex_ridge_list[i].second;
    D_(2, i) = ridge.x() / mass;
    D_(5, i) = ridge.y() / mass;
  }
}

LinearMpcXY::LinearMpcXY(double mass, double horizon_dt, QpSolverCollection::QpSolverType qp_solver_type)
: mass_(mass), horizon_dt_(horizon_dt), force_range_(3.0, 3.0 * mass * constants::g)
{
  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

Eigen::VectorXd LinearMpcXY::planOnce(const std::function<MotionParam(double)> & motion_param_func,
                                      const std::function<RefData(double)> & ref_data_func,
                                      const InitialParam & initial_param,
                                      const std::pair<double, double> & horizon_time_range,
                                      const WeightParam & weight_param)
{
  // Set model_list and ref_output_seq
  int horizon_size = static_cast<int>((horizon_time_range.second - horizon_time_range.first) / horizon_dt_);
  std::vector<std::shared_ptr<_StateSpaceModel>> model_list(horizon_size);
  Eigen::VectorXd ref_output_seq(horizon_size * RefData::outputDim());
  for(int i = 0; i < horizon_size; i++)
  {
    double t = horizon_time_range.first + i * horizon_dt_;
    model_list[i] = std::make_shared<Model>(mass_, motion_param_func(t));
    model_list[i]->calcDiscMatrix(horizon_dt_);
    ref_output_seq.segment<RefData::outputDim()>(i * RefData::outputDim()) = ref_data_func(t).toOutput(mass_);
  }

  // Calculate optimal force
  return procOnce(model_list, initial_param.toState(mass_), ref_output_seq, weight_param);
}

void LinearMpcXY::planLoop(const std::function<MotionParam(double)> & motion_param_func,
                           const std::function<RefData(double)> & ref_data_func,
                           const InitialParam & initial_param,
                           const std::pair<double, double> & motion_time_range,
                           double horizon_duration,
                           double sim_dt,
                           const WeightParam & weight_param)
{
  int seq_len = static_cast<int>((motion_time_range.second - motion_time_range.first) / sim_dt);
  int horizon_size = static_cast<int>(horizon_duration / horizon_dt_);

  // Loop
  double current_t = motion_time_range.first;
  StateDimVector current_x = initial_param.toState(mass_);
  motion_data_seq_.resize(seq_len);
  for(int i = 0; i < seq_len; i++)
  {
    // Set model_list and ref_output_seq
    std::vector<std::shared_ptr<_StateSpaceModel>> model_list(horizon_size);
    Eigen::VectorXd ref_output_seq(horizon_size * RefData::outputDim());
    const auto & current_motion_param = motion_param_func(current_t);
    const auto & current_ref_data = ref_data_func(current_t);
    for(int i = 0; i < horizon_size; i++)
    {
      double t = current_t + i * horizon_dt_;
      const auto & motion_param = (i == 0 ? current_motion_param : motion_param_func(t));
      const auto & ref_data = (i == 0 ? current_ref_data : ref_data_func(t));

      model_list[i] = std::make_shared<Model>(mass_, motion_param);
      model_list[i]->calcDiscMatrix(horizon_dt_);

      ref_output_seq.segment<RefData::outputDim()>(i * RefData::outputDim()) = ref_data.toOutput(mass_);
    }
    const auto & current_model = model_list[0];

    // Set sim_model
    const auto & sim_model = std::make_shared<SimModel>(mass_, current_motion_param);
    sim_model->calcDiscMatrix(sim_dt);

    // Calculate optimal force
    const Eigen::VectorXd & opt_force_seq = procOnce(model_list, current_x, ref_output_seq, weight_param);

    // Save current data
    Eigen::VectorXd current_u(current_model->inputDim());
    current_u << opt_force_seq.head(current_model->inputDim());
    const auto & current_sim_output = sim_model_->observEq(current_x, current_u);
    auto & current_motion_data = motion_data_seq_[i];
    current_motion_data.time = current_t;
    current_motion_data.ref_pos = current_ref_data.pos;
    current_motion_data.ref_vel = current_ref_data.vel;
    current_motion_data.ref_acc << 0.0, 0.0;
    current_motion_data.ref_angular_momentum = current_ref_data.angular_momentum;
    current_motion_data.planned_pos << current_sim_output[0], current_sim_output[3];
    current_motion_data.planned_vel << current_sim_output[1], current_sim_output[4];
    current_motion_data.planned_acc << current_sim_output[2], current_sim_output[5];
    current_motion_data.planned_force = current_u;
    current_motion_data.planned_angular_momentum = current_sim_output.tail<2>();

    // Simulate one step
    current_t += sim_dt;
    current_x = sim_model_->stateEqDisc(current_x, current_u);
  }
}

Eigen::VectorXd LinearMpcXY::procOnce(const std::vector<std::shared_ptr<_StateSpaceModel>> & model_list,
                                      const StateDimVector & current_x,
                                      const Eigen::VectorXd & ref_output_seq,
                                      const WeightParam & weight_param)
{
  // Calculate sequential extension
  int horizon_size = model_list.size();
  VariantSequentialExtension<state_dim_> seq_ext(model_list, false);

  // Set QP objective coefficients
  int total_input_dim = seq_ext.totalInputDim();
  if(!(qp_coeff_.dim_var_ == total_input_dim && qp_coeff_.dim_eq_ == horizon_size && qp_coeff_.dim_ineq_ == 0))
  {
    qp_coeff_.setup(total_input_dim, horizon_size, 0);
  }
  const Eigen::VectorXd & output_weight = weight_param.outputWeight(model_list.size());
  qp_coeff_.obj_mat_.noalias() = seq_ext.B_seq_.transpose() * output_weight.asDiagonal() * seq_ext.B_seq_;
  // Diagonal matrix cannot be added at the same time.
  // See https://forum.kde.org/viewtopic.php?f=74&t=136617#p365547
  qp_coeff_.obj_mat_ += weight_param.inputWeight(total_input_dim).asDiagonal();
  qp_coeff_.obj_vec_.noalias() = -1 * seq_ext.B_seq_.transpose() * output_weight.asDiagonal()
                                 * (ref_output_seq - seq_ext.A_seq_ * current_x - seq_ext.E_seq_);

  // Set QP constraint coefficients
  int accum_input_dim = 0;
  for(int i = 0; i < horizon_size; i++)
  {
    const auto & model = std::dynamic_pointer_cast<Model>(model_list[i]);

    for(size_t j = 0; j < model->motion_param_.vertex_ridge_list.size(); j++)
    {
      const auto & ridge = model->motion_param_.vertex_ridge_list[j].second;
      qp_coeff_.eq_mat_(i, accum_input_dim + j) = ridge.z();
    }
    qp_coeff_.eq_vec_[i] = model->motion_param_.total_force_z;

    accum_input_dim += model->inputDim();
  }
  qp_coeff_.x_min_.setConstant(force_range_.first);
  qp_coeff_.x_max_.setConstant(force_range_.second);

  // Solve QP
  return qp_solver_->solve(qp_coeff_);
}

void LinearMpcXY::dumpMotionDataSeq(const std::string & file_path, bool print_command) const
{
  std::ofstream ofs(file_path);
  ofs << "time ref_pos_x ref_pos_y ref_vel_x ref_vel_y planned_pos_x planned_pos_y planned_vel_x planned_vel_y "
         "planned_acc_x planned_acc_y planned_force"
      << std::endl;
  for(const auto & motion_data : motion_data_seq_)
  {
    motion_data.dump(ofs);
  }
  if(print_command)
  {
    std::cout << "Run the following commands in gnuplot:\n"
              << "  set key autotitle columnhead\n"
              << "  set key noenhanced\n"
              << "  plot \"" << file_path << "\" u 1:2 w lp, \"\" u 1:6 w lp\n"
              << "  replot \"" << file_path << "\" u 1:3 w lp, \"\" u 1:7 w lp\n";
  }
}
