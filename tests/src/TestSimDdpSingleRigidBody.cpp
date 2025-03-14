/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <ros/ros.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_srvs/Empty.h>

#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>

#include <CCC/Constants.h>
#include <CCC/DdpSingleRigidBody.h>

#include "ContactManager.h"
#include "SimModels.h"

class TestSimDdpSingleRigidBody
{
public:
  TestSimDdpSingleRigidBody()
  {
    // Setup ROS
    control_pub_ = nh_.advertise<std_msgs::Float64MultiArray>("control", 1);
    state_sub_ = nh_.subscribe("state", 1, &TestSimDdpSingleRigidBody::stateCallback, this);
    forward_srv_ = nh_.advertiseService("/forward", &TestSimDdpSingleRigidBody::forwardCallback, this);
    jump_srv_ = nh_.advertiseService("/jump", &TestSimDdpSingleRigidBody::jumpCallback, this);
    tilt_srv_ = nh_.advertiseService("/tilt", &TestSimDdpSingleRigidBody::tiltCallback, this);
  }

  void run()
  {
    double horizon_dt = 0.03; // [sec]
    double horizon_duration = 3.0; // [sec]
    int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);

    // Setup DDP
    CCC::DdpSingleRigidBody::WeightParam weight_param;
    weight_param.running_pos << 1.0, 1.0, 10.0;
    weight_param.running_ori << 0.5, 0.5, 0.5;
    weight_param.terminal_pos << 1.0, 1.0, 10.0;
    weight_param.terminal_ori << 0.5, 0.5, 0.5;
    CCC::DdpSingleRigidBody ddp(mass_, horizon_dt, horizon_steps, weight_param);
    ddp.ddp_solver_->config().max_iter = 1;
    initial_param_.pos = Eigen::Vector3d(0.0, 0.0, 1.0);

    // Setup contact
    std::function<CCC::DdpSingleRigidBody::MotionParam(double)> motion_param_func = [this](double t)
    {
      CCC::DdpSingleRigidBody::MotionParam motion_param;
      Eigen::Vector2d contact_pos = Eigen::Vector2d::Zero();
      if(forward_duration_ && (*forward_duration_)[0] <= t && t <= (*forward_duration_)[1])
      {
        contact_pos.x() += forward_dist_;
      }
      if(!(jump_duration_ && (*jump_duration_)[0] <= t && t <= (*jump_duration_)[1]))
      {
        motion_param.contact_list.push_back(
            makeContactFromRect({contact_pos + Eigen::Vector2d(-0.5, -0.5), contact_pos + Eigen::Vector2d(0.5, 0.5)}));
      }
      motion_param.inertia_mat.diagonal() = moment_of_inertia_;
      return motion_param;
    };
    std::function<CCC::DdpSingleRigidBody::RefData(double)> ref_data_func = [this](double t)
    {
      CCC::DdpSingleRigidBody::RefData ref_data;
      ref_data.pos << 0.0, 0.0, 1.0;
      if(forward_duration_ && (*forward_duration_)[0] <= t && t <= (*forward_duration_)[1])
      {
        ref_data.pos.x() += forward_dist_;
      }
      if(jump_duration_ && (*jump_duration_)[0] <= t && t <= (*jump_duration_)[1])
      {
        ref_data.pos.z() += jump_height_;
      }
      ref_data.ori.setZero();
      if(tilt_duration_ && (*tilt_duration_)[0] <= t && t <= (*tilt_duration_)[1])
      {
        ref_data.ori.x() += tilt_angle_;
        ref_data.ori.y() += tilt_angle_;
      }
      return ref_data;
    };

    // Run control loop
    ros::Rate rate(200);
    while(ros::ok())
    {
      t_ = ros::Time::now().toSec();

      if(forward_duration_ && (*forward_duration_)[1] < t_)
      {
        forward_duration_.reset();
      }
      if(jump_duration_ && (*jump_duration_)[1] < t_)
      {
        jump_duration_.reset();
      }
      if(tilt_duration_ && (*tilt_duration_)[1] < t_)
      {
        tilt_duration_.reset();
      }

      ros::spinOnce();

      // Plan
      if(!initial_param_.u_list.empty())
      {
        for(int i = 0; i < ddp.ddp_solver_->config().horizon_steps; i++)
        {
          double tmp_time = t_ + i * ddp.ddp_problem_->dt();
          int input_dim = ddp.ddp_problem_->inputDim(tmp_time);
          if(initial_param_.u_list[i].size() != input_dim)
          {
            initial_param_.u_list[i].setZero(input_dim);
          }
        }
      }
      Eigen::VectorXd planned_force_scales = ddp.planOnce(motion_param_func, ref_data_func, initial_param_, t_);

      // Publish
      std_msgs::Float64MultiArray msg;
      const auto & motion_param = motion_param_func(t_);
      for(const auto & contact : motion_param.contact_list)
      {
        int wrenchRatioIdx = 0;
        for(const auto & vertexWithRidge : contact->vertexWithRidgeList_)
        {
          const Eigen::Vector3d & vertex = vertexWithRidge.vertex;
          const std::vector<Eigen::Vector3d> & ridgeList = vertexWithRidge.ridgeList;

          Eigen::Vector3d vertexForce = Eigen::Vector3d::Zero();
          for(const auto & ridge : ridgeList)
          {
            vertexForce += planned_force_scales(wrenchRatioIdx) * ridge;
            wrenchRatioIdx++;
          }

          std::vector<double> vertexVec(vertex.data(), vertex.data() + 3);
          std::vector<double> vertexForceVec(vertexForce.data(), vertexForce.data() + 3);
          msg.data.insert(msg.data.end(), vertexVec.begin(), vertexVec.end());
          msg.data.insert(msg.data.end(), vertexForceVec.begin(), vertexForceVec.end());
        }
      }
      control_pub_.publish(msg);

      rate.sleep();
    }
  }

protected:
  void stateCallback(const std_msgs::Float64MultiArray::ConstPtr & msg)
  {
    const CCC::DdpSingleRigidBody::DdpProblem::StateDimVector & state =
        Eigen::Map<const CCC::DdpSingleRigidBody::DdpProblem::StateDimVector>(msg->data.data());
    initial_param_.pos = state.segment<3>(0);
    initial_param_.ori = state.segment<3>(3);
    initial_param_.linear_vel = state.segment<3>(6);
    initial_param_.angular_vel = state.segment<3>(9);
  }

  bool forwardCallback(std_srvs::Empty::Request &, // req
                       std_srvs::Empty::Response & // res
  )
  {
    if(forward_duration_)
    {
      return false;
    }

    forward_duration_ = std::make_shared<std::array<double, 2>>();
    (*forward_duration_)[0] = t_ + 2.0;
    (*forward_duration_)[1] = t_ + 4.0;

    return true;
  }

  bool jumpCallback(std_srvs::Empty::Request &, // req
                    std_srvs::Empty::Response & // res
  )
  {
    if(jump_duration_)
    {
      return false;
    }

    jump_duration_ = std::make_shared<std::array<double, 2>>();
    (*jump_duration_)[0] = t_ + 2.0;
    (*jump_duration_)[1] = t_ + 2.4;

    return true;
  }

  bool tiltCallback(std_srvs::Empty::Request &, // req
                    std_srvs::Empty::Response & // res
  )
  {
    if(tilt_duration_)
    {
      return false;
    }

    tilt_duration_ = std::make_shared<std::array<double, 2>>();
    (*tilt_duration_)[0] = t_ + 2.0;
    (*tilt_duration_)[1] = t_ + 5.0;

    return true;
  }

protected:
  double t_ = 0.0;

  double mass_ = 10.0; // [kg]
  Eigen::Vector3d moment_of_inertia_ =
      Eigen::Vector3d(0.6166666666666667, 0.48333333333333334, 0.2833333333333333); // [kg m^2]

  CCC::DdpSingleRigidBody::InitialParam initial_param_;

  const double forward_dist_ = 1.5; // [m]
  const double jump_height_ = 0.4; // [m]
  const double tilt_angle_ = 1.0; // [rad]
  std::shared_ptr<std::array<double, 2>> forward_duration_;
  std::shared_ptr<std::array<double, 2>> jump_duration_;
  std::shared_ptr<std::array<double, 2>> tilt_duration_;

  ros::NodeHandle nh_;
  ros::Publisher control_pub_;
  ros::Subscriber state_sub_;
  ros::ServiceServer forward_srv_;
  ros::ServiceServer jump_srv_;
  ros::ServiceServer tilt_srv_;
};

TEST(TestSimDdpSingleRigidBody, Test1)
{
  TestSimDdpSingleRigidBody test;

  test.run();
}

int main(int argc, char ** argv)
{
  // Setup ROS
  ros::init(argc, argv, "test_sim_ddp_srb");

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
