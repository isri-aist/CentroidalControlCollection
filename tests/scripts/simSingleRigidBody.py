#! /usr/bin/env python3

import numpy as np
import pybullet
import pybullet_data
import unittest
import sys

import rclpy
from tf_transformations import quaternion_from_euler, euler_from_quaternion
from std_msgs.msg import Float64MultiArray
import launch
import launch_testing
import pytest
from launch_ros.actions import Node

@pytest.mark.launch_test
def generate_test_description():
    test_sim_ddp_srb = Node(
        package='centroidal_control_collection',
        executable='TestSimDdpSingleRigidBody',
        name='test_sim_ddp_srb',
        output='screen'
    )

    content = {}

    return (
        launch.LaunchDescription([
            test_sim_ddp_srb,
            launch_testing.actions.ReadyToTest(),
        ]),
        content
    )

class SimSingleRigidBody(object):
    def __init__(self, enable_gui=True):
        self.node = rclpy.create_node("sim")
        # Instantiate simulator
        if enable_gui:
            pybullet.connect(pybullet.GUI)
        else:
            pybullet.connect(pybullet.DIRECT)

        # Set simulation parameters
        self.dt = 0.002 # [sec]
        pybullet.setTimeStep(self.dt)
        pybullet.setGravity(0, 0, -9.80665) # [m/s^2]

        # Set debug parameters
        pybullet.configureDebugVisualizer(pybullet.COV_ENABLE_GUI, 0)

        # Setup models
        pybullet.setAdditionalSearchPath(pybullet_data.getDataPath())

        ## Setup floor
        pybullet.loadURDF("plane100.urdf")

        ## Setup body
        self.box_half_scale = np.array([0.15, 0.25, 0.35]) # [m]
        box_col_shape_idx = pybullet.createCollisionShape(pybullet.GEOM_BOX,
                                                          halfExtents=self.box_half_scale)
        box_mass = 10.0 # [kg]
        self.body_uid = pybullet.createMultiBody(baseMass=box_mass,
                                                 baseCollisionShapeIndex=box_col_shape_idx,
                                                 baseVisualShapeIndex=-1,
                                                 basePosition=[0.0, 0.0, 1.0],
                                                 baseOrientation=[0.0, 0.0, 0.0, 1.0],
                                                 baseInertialFramePosition=[0.0, 0.0, 0.0],
                                                 baseInertialFrameOrientation=[0.0, 0.0, 0.0, 1.0])
        pybullet.changeVisualShape(objectUniqueId=self.body_uid,
                                   linkIndex=-1,
                                   rgbaColor=[0.0, 1.0, 0.0, 0.8])
        info = pybullet.getDynamicsInfo(bodyUniqueId=self.body_uid,
                                        linkIndex=-1)
        print(info)

        # Setup variables
        self.force_scale = 0.01
        self.pos_force_list = []
        self.force_line_uid_list = []

        # Setup ROS
        self.state_pub = self.node.create_publisher(Float64MultiArray, "state", 1)
        self.control_sub = self.node.create_subscription(Float64MultiArray, "control", self.callback, 1)

    def runOnce(self):
        """"Run simulation step once."""
        # Apply force
        c, _ = pybullet.getBasePositionAndOrientation(bodyUniqueId=self.body_uid)
        wrench = np.zeros(6)
        for pos, force in self.pos_force_list:
            wrench += np.array([force, np.cross(pos - c, force)]).flatten()
        pybullet.applyExternalForce(objectUniqueId=self.body_uid,
                                    linkIndex=-1,
                                    forceObj=wrench[0:3],
                                    posObj=c,
                                    flags=pybullet.WORLD_FRAME)
        pybullet.applyExternalTorque(objectUniqueId=self.body_uid,
                                     linkIndex=-1,
                                     torqueObj=wrench[3:6],
                                     flags=pybullet.WORLD_FRAME)

        # Visualize force line
        force_idx = 0
        for pos, force in self.pos_force_list:
            if force_idx >= len(self.force_line_uid_list):
                self.force_line_uid_list.append(-1)
            self.force_line_uid_list[force_idx] = \
                pybullet.addUserDebugLine(lineFromXYZ=pos,
                                          lineToXYZ=pos + self.force_scale * force,
                                          lineColorRGB=[1, 0, 0],
                                          lineWidth=5.0,
                                          replaceItemUniqueId=self.force_line_uid_list[force_idx])
            force_idx += 1
        # TODO: The following code should be necessary to remove unnecessary visible lines,
        #       but comment it out because it slows down the simulation
        # if force_idx < len(self.force_line_uid_list):
        #     for i in range(force_idx, len(self.force_line_uid_list)):
        #         pybullet.removeUserDebugItem(self.force_line_uid_list[i])
        #     self.force_line_uid_list = self.force_line_uid_list[:len(self.pos_force_list)]

        # Publish state
        state_msg = Float64MultiArray()
        state_msg.data = self.getState()
        self.state_pub.publish(state_msg)

        # Process simulation step
        pybullet.stepSimulation()

    def getState(self):
        """"Get state [c, alpha, v, omega]."""
        c, quat = pybullet.getBasePositionAndOrientation(bodyUniqueId=self.body_uid)
        alpha = euler_from_quaternion(quat, axes="rzyx")
        v, omega = pybullet.getBaseVelocity(bodyUniqueId=self.body_uid)
        return np.array([c, alpha, v, omega]).flatten().tolist()

    def setState(self, state):
        """Set state [c, alpha, v, omega]."""
        pybullet.resetBasePositionAndOrientation(bodyUniqueId=self.body_uid,
                                                 posObj=state[0:3],
                                                 ornObj=quaternion_from_euler(*state[3:6], axes="rzyx"))
        pybullet.resetBaseVelocity(objectUniqueId=self.body_uid,
                                   linearVelocity=state[6:9],
                                   angularVelocity=state[9:12])

    def callback(self, msg):
        contact_num = int(len(msg.data) / 6)
        self.pos_force_list = []
        for i in range(contact_num):
            self.pos_force_list.append((np.array(msg.data[6*i+0:6*i+3]), np.array(msg.data[6*i+3:6*i+6])))

def demo():
    sim = SimSingleRigidBody(True)

    t = 0.0 # [sec]
    rate = sim.node.create_rate(1.0 / sim.dt)
    while pybullet.isConnected():
        # Run simulation step
        sim.runOnce()

        # Sleep and increment time
        rate.sleep()
        t += sim.dt

class TestSimSingleRigidBody(unittest.TestCase):
    def __init__(self, *args):
        super().__init__(*args)

    def test(self):
        rclpy.init(args=sys.argv)
        demo()