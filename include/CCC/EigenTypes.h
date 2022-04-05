#pragma once

#include <Eigen/Core>

namespace Eigen
{
template<typename T>
using Vector1 = Matrix<T, 1, 1>;

template<typename T>
using Vector6 = Matrix<T, 6, 1>;

typedef Vector1<double> Vector1d;
typedef Vector6<double> Vector6d;
} // namespace Eigen
