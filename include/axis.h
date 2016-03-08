#ifndef UNIT_AXES_H
#define UNIT_AXES_H

#include <cmath>
#include <Eigen/Dense>

namespace axes
{

const double EPS = 1e-9;

///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis2D
///////////////////////////////////////////////////////////////////////////////////////////////////

class UnitAxis2D
{
public:
  // Constructors
  UnitAxis2D(double, double);
  UnitAxis2D(double phi) : UnitAxis2D(std::cos(phi), std::sin(phi)) {}
  UnitAxis2D() : UnitAxis2D(1.0, 0.0) {}
  UnitAxis2D(const Eigen::Vector2d& v) : UnitAxis2D(v(0), v(1)) {}
  UnitAxis2D(const UnitAxis2D& orig) : lambda_(orig.lambda_), kappa_(orig.kappa_) {}
  // Introspection
  double          angle() const;
  double          kappa() const;
  double          lambda() const;
  Eigen::Vector2d vector() const;
  // Mathematical methods
  UnitAxis2D      inv() const;
  double          log() const;
  // Operators
  double          operator[](std::size_t) const;
  Eigen::Matrix2d operator+() const;
  Eigen::Matrix2d operator-() const;
  UnitAxis2D&     operator=(const UnitAxis2D&);
private:
  double          lambda_;
  double          kappa_;
};

// Functions acting on UnitAxis2D
std::ostream &operator<<(std::ostream&, const UnitAxis2D&);
bool         operator==(const UnitAxis2D&, const UnitAxis2D&);
UnitAxis2D   operator+(const UnitAxis2D&, const UnitAxis2D&);
double       boxminus(const UnitAxis2D&, const UnitAxis2D&);
UnitAxis2D   boxplus(const UnitAxis2D&, const Eigen::Vector2d&);
double       distance(const UnitAxis2D&, const UnitAxis2D&);
double       dot(const UnitAxis2D&, const UnitAxis2D&);

///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis3D
///////////////////////////////////////////////////////////////////////////////////////////////////

class UnitAxis3D
{
public:
  // Constructors
  UnitAxis3D(double, double,  double);
  UnitAxis3D() : UnitAxis3D(1.0, 0.0, 0.0) {}
  UnitAxis3D(const Eigen::Vector2d&);
  UnitAxis3D(const Eigen::Vector3d& v) : UnitAxis3D(v(0), v(1), v(2)) {}
  UnitAxis3D(double lambda, const Eigen::Vector2d& kappa) : UnitAxis3D(lambda, kappa(0), kappa(1)) {}
  UnitAxis3D(const UnitAxis3D& orig) : lambda_(orig.lambda_), kappa_(orig.kappa_) {}
  // Instrospection
  Eigen::Vector2d direction() const;
  double          inclination() const;
  Eigen::Vector2d kappa() const;
  double          lambda() const;
  Eigen::Vector3d vector() const;
   // Mathematical methods
  UnitAxis3D      inv() const;
  Eigen::Vector2d log() const;
  // Operators
  double          operator[](std::size_t) const;
  Eigen::Matrix3d operator+() const;
  Eigen::Matrix3d operator-() const;
  UnitAxis3D&     operator=(const UnitAxis3D&);
private:
  double          lambda_;
  Eigen::Vector2d kappa_;
};

// Functions acting on UnitAxis3D
std::ostream    &operator<<(std::ostream&, const UnitAxis3D&);
bool            operator==(const UnitAxis3D&, const UnitAxis3D&);
UnitAxis3D      operator+(const UnitAxis3D&, const UnitAxis3D&);
Eigen::Vector2d boxminus(const UnitAxis3D&, const UnitAxis3D&);
UnitAxis3D      boxplus(const UnitAxis3D&, const Eigen::Vector3d&);
double          distance(const UnitAxis3D&, const UnitAxis3D&);
double          dot(const UnitAxis3D&, const UnitAxis3D&);
UnitAxis3D      exp(const Eigen::Vector2d&);

} // namespace axes

#endif
