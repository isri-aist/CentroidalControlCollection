#pragma once

#include <Eigen/Core>

namespace Eigen
{
template<typename T>
using Vector1 = Matrix<T, 1, 1>;

typedef Vector1<double> Vector1d;
}
