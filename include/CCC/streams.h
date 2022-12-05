#pragma once

#ifdef CCC_STANDALONE
#  include <iostream>
#  define CCC_ERROR_STREAM(x) std::cerr << x << "\n"
#  define CCC_WARN_STREAM(x) std::cerr << x << "\n"
#  define CCC_INFO_STREAM(x) std::cerr << x << "\n"
#else
#  include <ros/console.h>
#  define CCC_ERROR_STREAM ROS_ERROR_STREAM
#  define CCC_WARN_STREAM ROS_WARN_STREAM
#  define CCC_INFO_STREAM ROS_INFO_STREAM
#endif
