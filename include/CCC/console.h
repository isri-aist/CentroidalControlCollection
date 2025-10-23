#pragma once

#ifdef CCC_STANDALONE
#  include <iostream>
#  define CCC_ERROR_STREAM(x) std::cerr << x << "\n"
#  define CCC_WARN_STREAM(x) std::cerr << x << "\n"
#  define CCC_INFO_STREAM(x) std::cout << x << "\n"
#else
#  include <rclcpp/rclcpp.hpp>
#  define CCC_ERROR_STREAM(msg) RCLCPP_ERROR_STREAM(rclcpp::get_logger("CentroidalControlCollection"), msg)
#  define CCC_WARN_STREAM(msg) RCLCPP_WARN_STREAM(rclcpp::get_logger("CentroidalControlCollection"), msg)
#  define CCC_INFO_STREAM(msg) RCLCPP_INFO_STREAM(rclcpp::get_logger("CentroidalControlCollection"), msg)
#endif
