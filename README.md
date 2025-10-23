**This is the branch for ROS1.**

# [CentroidalControlCollection](https://github.com/isri-aist/CentroidalControlCollection)
Collection of centroidal control for legged robots

[![CI-standalone](https://github.com/isri-aist/CentroidalControlCollection/actions/workflows/ci-standalone.yaml/badge.svg)](https://github.com/isri-aist/CentroidalControlCollection/actions/workflows/ci-standalone.yaml)
[![CI-catkin](https://github.com/isri-aist/CentroidalControlCollection/actions/workflows/ci-catkin.yaml/badge.svg)](https://github.com/isri-aist/CentroidalControlCollection/actions/workflows/ci-catkin.yaml)
[![Documentation](https://img.shields.io/badge/doxygen-online-brightgreen?logo=read-the-docs&style=flat)](https://isri-aist.github.io/CentroidalControlCollection/)

## Install

### Requirements
- Compiler supporting C++17
- Tested on `Ubuntu 20.04 / ROS Noetic` and `Ubuntu 18.04 / ROS Melodic`

### Dependencies
This package depends on
- [mc_rtc](https://jrl-umi3218.github.io/mc_rtc)

This package also depends on the following packages. However, manual installation is unnecessary when this package is installed using `wstool` as described in [Installation procedure](#installation-procedure).
- [QpSolverCollection](https://github.com/isri-aist/QpSolverCollection)
- [ForceControlCollection](https://github.com/isri-aist/ForceControlCollection)
- [NMPC](https://github.com/isri-aist/NMPC)

### Preparation
1. (Skip if ROS is already installed.) Install ROS. See [here](http://wiki.ros.org/ROS/Installation) for details.
```bash
$ export ROS_DISTRO=melodic
$ sudo sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list'
$ wget http://packages.ros.org/ros.key -O - | sudo apt-key add -
$ sudo apt-get update
$ sudo apt-get install ros-${ROS_DISTRO}-ros-base python-catkin-tools python-rosdep
```

2. (Skip if mc_rtc is already installed.) Install mc_rtc. See [here](https://jrl-umi3218.github.io/mc_rtc/tutorials/introduction/installation-guide.html) for details.
```bash
$ curl -1sLf 'https://dl.cloudsmith.io/public/mc-rtc/stable/setup.deb.sh' | sudo -E bash
$ sudo apt-get install libmc-rtc-dev mc-rtc-utils ros-${ROS_DISTRO}-mc-rtc-plugin ros-${ROS_DISTRO}-mc-rtc-rviz-panel libeigen-qld-dev
```

### Installation procedure
1. Setup catkin workspace.
```bash
$ mkdir -p ~/ros/ws_ccc/src
$ cd ~/ros/ws_ccc
$ wstool init src
$ wstool set -t src isri-aist/QpSolverCollection git@github.com:isri-aist/QpSolverCollection.git --git -y
$ wstool set -t src isri-aist/ForceControlCollection git@github.com:isri-aist/ForceControlCollection.git --git -y
$ wstool set -t src isri-aist/NMPC git@github.com:isri-aist/NMPC.git --git -y
$ wstool set -t src isri-aist/CentroidalControlCollection git@github.com:isri-aist/CentroidalControlCollection.git --git -y
$ wstool update -t src
```

2. Install dependent packages.
```bash
$ source /opt/ros/${ROS_DISTRO}/setup.bash
$ rosdep install -y -r --from-paths src --ignore-src
```

3. Build a package.
```bash
$ catkin build centroidal_control_collection -DCMAKE_BUILD_TYPE=RelWithDebInfo --catkin-make-args all tests
```

## Examples
Make sure that it is built with `--catkin-make-args tests` option.

### Methods based on bipedal dynamics
The CoM and ZMP trajectories are planned according to the ZMP reference trajectory and the ZMP region boundaries as inputs, which are determined from a given footstep sequence (i.e., the position and timing of the foot landings). The CoM velocity is jumped by emulating a disturbance during motion.

#### [PreviewControlZmp](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1PreviewControlZmp.html)
- Shuuji Kajita, et al. Biped walking pattern generation by using preview control of zero-moment point. ICRA, 2003.

```bash
$ rosrun centroidal_control_collection TestPreviewControlZmp
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method PreviewControlZmp
```

![PreviewControlZmp](doc/images/PreviewControlZmp.png)

#### [DdpZmp](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1DdpZmp.html)
- S Feng, et al. Optimization‐based full body control for the darpa robotics challenge. Journal of field robotics, 2015.

```bash
$ rosrun centroidal_control_collection TestDdpZmp
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method DdpZmp
```

![DdpZmp](doc/images/DdpZmp.png)

#### [DcmTracking](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1DcmTracking.html)
- J Englsberger, et al. Three-dimensional bipedal walking control using divergent component of motion. IROS, 2013.

```bash
$ rosrun centroidal_control_collection TestDcmTracking
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method DcmTracking
```

![DcmTracking](doc/images/DcmTracking.png)

#### [FootGuidedControl](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1FootGuidedControl.html)
- T Sugihara, et al. Foot-guided agile control of a biped robot through ZMP manipulation. IROS, 2017.
- Y Kojio, et al. Unified balance control for biped robots including modification of footsteps with angular momentum and falling detection based on capturability. IROS, 2019.

```bash
$ rosrun centroidal_control_collection TestFootGuidedControl
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method FootGuidedControl
```

![FootGuidedControl](doc/images/FootGuidedControl.png)

#### [StepMpc](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1StepMpc.html)
- S Xin, et al. Online relative footstep optimization for legged robots dynamic walking using discrete-time model predictive control. IROS, 2019.

```bash
$ rosrun centroidal_control_collection TestStepMpc
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method StepMpc
```

![StepMpc](doc/images/StepMpc.png)

#### [LinearMpcZmp](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1LinearMpcZmp.html)
- PB Wieber. Trajectory Free Linear Model Predictive Control for Stable Walking in the Presence of Strong Perturbations. Humanoids, 2006.

```bash
$ rosrun centroidal_control_collection TestLinearMpcZmp
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method LinearMpcZmp
```

![LinearMpcZmp](doc/images/LinearMpcZmp.png)

#### [IntrinsicallyStableMpc](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1IntrinsicallyStableMpc.html)
- N Scianca, et al. Intrinsically Stable MPC for Humanoid Gait Generation. Humanoids, 2016.

```bash
$ rosrun centroidal_control_collection TestIntrinsicallyStableMpc
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method IntrinsicallyStableMpc
```

![IntrinsicallyStableMpc](doc/images/IntrinsicallyStableMpc.png)

#### [SingularPreviewControlZmp](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1SingularPreviewControlZmp.html)
- J Urata, et al. Online Decision of Foot Placement using Singular LQ Preview Regulation. Humanoids, 2011.

```bash
$ rosrun centroidal_control_collection TestSingularPreviewControlZmp
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method SingularPreviewControlZmp
```

![SingularPreviewControlZmp](doc/images/SingularPreviewControlZmp.png)

#### Plotting all methods

```bash
$ roscd centroidal_control_collection
$ catkin bt --no-deps --catkin-make-args run_tests
$ rosrun centroidal_control_collection plotTestZmpBasedMethodResults.py --method All --plot-comp-time
```

You will get a plot like [this one](https://www.dropbox.com/s/t2h9qu5tzfon4mq/plotTestZmpBasedMethodResultsAll.pdf?dl=0) that shows all the methods in one sheet.

### Methods based on centroidal dynamics
Centroidal trajectories (i.e., CoM and linear/angular momentum trajectories) are planned from the contact sequence.

#### [LinearMpcZ](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1LinearMpcZ.html)

```bash
$ rosrun centroidal_control_collection TestLinearMpcZ
```

#### [LinearMpcXY](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1LinearMpcXY.html)
- H Audren, et al. Model preview control in multi-contact motion-application to a humanoid robot. IROS, 2014.
- 長阪憲一郎, et al. 接触拘束を考慮可能なマルチコンタクト対応スタビライザと一般化逆動力学による人型ロボットの全身制御. ロボティクスシンポジア予稿集, 2012.

```bash
$ rosrun centroidal_control_collection TestLinearMpcXY
```

#### [PreviewControlCentroidal](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1PreviewControlCentroidal.html)
- M Murooka, et al. Centroidal trajectory generation and stabilization based on preview control for humanoid multi-contact motion. RA-Letters, 2022.

```bash
$ rosrun centroidal_control_collection TestPreviewControlCentroidal
```

#### [DdpCentroidal](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1DdpCentroidal.html)

```bash
$ rosrun centroidal_control_collection TestDdpCentroidal --gtest_filter=*.PlanOnce
```

#### [DdpSingleRigidBody](https://isri-aist.github.io/CentroidalControlCollection/doxygen/classCCC_1_1DdpSingleRigidBody.html)

```bash
$ rosrun centroidal_control_collection TestDdpSingleRigidBody --gtest_filter=*.PlanOnce
```

If you catkin build with `-DENABLE_PYBULLET_TEST=ON` option, you can run the test by simulation on the GUI. `pybullet` is required to be installed.
```bash
# In terminal 1
$ rostest centroidal_control_collection TestSimDdpSingleRigidBody.test -t -r
# In terminal 2
$ rosrun rqt_service_caller rqt_service_caller
```

https://github.com/isri-aist/CentroidalControlCollection/assets/6636600/7f70f728-f76d-49b3-a885-d6e2ab5d8594

## Integration into controller
Some of the methods implemented in this library are available in the humanoid controller [BaselineWalkingController](https://github.com/isri-aist/BaselineWalkingController) and [MultiContactController](https://github.com/isri-aist/MultiContactController).
