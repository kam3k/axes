#include "axis.h"
#include <exception>

namespace axes
{
  
///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis2D
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructors

UnitAxis2D::UnitAxis2D(double lambda, double kappa)
{
  // Very small values are exactly zero
  lambda = std::fabs(lambda) < EPS ? 0.0 : lambda;
  kappa = std::fabs(kappa) < EPS ? 0.0 : kappa;
  // Set lambda = 1.0 if both are zero
  lambda = (std::fabs(lambda) < EPS && std::fabs(kappa) < EPS) ? 1.0 : lambda;
  // Change signs if lambda is negative (order here is important!)
  kappa = lambda < 0.0 ? -kappa : kappa;
  lambda = lambda < 0.0 ? -lambda : lambda;
  // Ensure if lambda = 0, kappa = 1 (note lambda would be exactly 0.0 at this point)
  kappa = lambda == 0.0 ? 1.0 : kappa;
  // Normalize parameters (if required)
  double square_sum = lambda * lambda + kappa * kappa;
  if (std::fabs(square_sum - 1.0) > EPS)
  {
    double normalizer = std::sqrt(square_sum);
    lambda_ = lambda / normalizer;
    kappa_ = kappa / normalizer;
  }
  else
  {
    lambda_ = lambda;
    kappa_ = kappa;
  }
}

// Introspection

double UnitAxis2D::angle() const
{
  return std::atan(kappa_ / lambda_);
}

double UnitAxis2D::kappa() const
{
  return kappa_;
}

double UnitAxis2D::lambda() const
{
  return lambda_;
}

Eigen::Vector2d UnitAxis2D::vector() const
{
  return Eigen::Vector2d(lambda_, kappa_);
}

// Mathematical Methods

UnitAxis2D UnitAxis2D::inv() const
{
  return UnitAxis2D(lambda_, -kappa_);
}

double UnitAxis2D::log() const
{
  return angle();
}

// Operators

double UnitAxis2D::operator[](std::size_t i) const
{
  if (i != 0 && i != 1)
  {
    throw std::out_of_range("UnitAxis2D::operator[]");
  }
  return i == 0 ? lambda_ : kappa_;
}

Eigen::Matrix2d UnitAxis2D::operator+() const
{
  Eigen::Matrix2d ret;
  ret << lambda_, -kappa_, kappa_, lambda_;
  return ret;
}

Eigen::Matrix2d UnitAxis2D::operator-() const
{
  Eigen::Matrix2d ret;
  ret << lambda_, kappa_, -kappa_, lambda_;
  return ret;
}

UnitAxis2D& UnitAxis2D::operator=(const UnitAxis2D& rhs)
{
  lambda_ = rhs.lambda();
  kappa_ = rhs.kappa();
  return *this;
}

// Functions acting on UnitAxis2D

std::ostream& operator<<(std::ostream& os, const UnitAxis2D& m)
{
  os << "lambda: " << m.lambda() << ", kappa: " << m.kappa();
  return os;
}

bool operator==(const UnitAxis2D& m, const UnitAxis2D& n)
{
  return std::fabs(m.lambda() - n.lambda()) < EPS && std::fabs(m.kappa() - n.kappa()) < EPS;
}

bool operator!=(const UnitAxis2D& m, const UnitAxis2D& n)
{
  return !(m == n);
}

UnitAxis2D operator+(const UnitAxis2D& m, const UnitAxis2D& n)
{
  return UnitAxis2D(+m * n.vector());
}

UnitAxis2D operator-(const UnitAxis2D& m, const UnitAxis2D& n)
{
  return UnitAxis2D(-m * n.vector());
}

double boxminus(const UnitAxis2D& m, const UnitAxis2D& n)
{
  return (n - m).log();
}

UnitAxis2D boxplus(const UnitAxis2D& m, double phi)
{
  return UnitAxis2D(m + UnitAxis2D(phi));
}

double distance(const UnitAxis2D& m, const UnitAxis2D& n)
{
  return std::acos(dot(m, n));
}

double dot(const UnitAxis2D& m, const UnitAxis2D& n)
{
  double x = std::fabs(m.lambda() * n.lambda() + m.kappa() * n.kappa());
  return x > 1.0 ? 1.0 : x; // numerical check to ensure x <= 1.0
}

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
