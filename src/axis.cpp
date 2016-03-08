#include "axis.h"

namespace axes
{

double UnitAxis2D::angle() const
{
  return 0.0;
}

Eigen::Vector2d UnitAxis3D::direction() const
{
  Eigen::Vector2d ret(0.0, 0.0);
  return ret;
}

double UnitAxis3D::inclination() const
{
  return 0.0;
}

} // namespace axes
