#include "axis.h"

namespace axes
{
  
///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis2D
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructors

UnitAxis2D::UnitAxis2D(double lambda, double kappa)
{
  // Both zero returns identity unit axis
  lambda = (std::fabs(lambda) < EPS && std::fabs(kappa) < EPS) ? 1.0 : lambda;
  // Very small values are exactly zero
  lambda = std::fabs(lambda) < EPS ? 0.0 : lambda;
  kappa = std::fabs(kappa) < EPS ? 0.0 : kappa;
  // Normalize parameters
  double normalizer = 1.0 / std::sqrt(lambda * lambda + kappa * kappa);
  lambda = lambda * normalizer;
  kappa = kappa * normalizer;
  // Ensure lambda >= 0, and kappa = 1.0 if lambda = 0
  lambda_ = lambda > 0.0 ? lambda : -lambda;
  kappa_ = lambda >= 0.0 ? kappa : -kappa;
}

// Introspection

double UnitAxis2D::angle() const
{
  return 0.0;
}

double UnitAxis2D::kappa() const
{
  return kappa_;
}

double UnitAxis2D::lambda() const
{
  return lambda_;
}

// Mathematical Methods

// Operators

///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis3D
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructors

// Introspection

Eigen::Vector2d UnitAxis3D::direction() const
{
  Eigen::Vector2d ret(0.0, 0.0);
  return ret;
}

double UnitAxis3D::inclination() const
{
  return 0.0;
}

// Mathematical Methods

// Operators

} // namespace axes
