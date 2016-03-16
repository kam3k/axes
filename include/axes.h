#ifndef UNIT_AXES_H
#define UNIT_AXES_H

#include <cmath>
#include <exception>
#include <Eigen/Dense>

namespace axes
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis2D
///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
class UnitAxis2D
{
public:
  // Typedefs
  typedef Eigen::Matrix<Scalar, 2, 1> Vector;
  typedef Eigen::Matrix<Scalar, 2, 2> Matrix;
  // Constructors
  UnitAxis2D(Scalar, Scalar);
  UnitAxis2D(Scalar phi) : UnitAxis2D(std::cos(phi), std::sin(phi)) {}
  UnitAxis2D() : UnitAxis2D(1, 0) {}
  UnitAxis2D(const Vector v) : UnitAxis2D(v(0), v(1)) {}
  UnitAxis2D(const UnitAxis2D& orig) : lambda_(orig.lambda_), kappa_(orig.kappa_) {}
  // Introspection
  /* Returns the angle of the UnitAxis */
  Scalar      angle() const {return std::atan(kappa_ / lambda_);}
  Scalar      kappa() const {return kappa_;}
  Scalar      lambda() const {return lambda_;}
  Vector      vector() const {return Vector(lambda_, kappa_);}
  // Mathematical methods
  UnitAxis2D  inv() const {return UnitAxis2D<Scalar>(lambda_, -kappa_);}
  Scalar      log() const {return angle();}
  // Operators
  Scalar      operator[](std::size_t) const;
  Matrix      operator+() const;
  Matrix      operator-() const;
  UnitAxis2D& operator=(const UnitAxis2D&);
private:
  Scalar      lambda_;
  Scalar      kappa_;
  Scalar      eps_ = std::numeric_limits<Scalar>::epsilon();
};

// Constructor

template <typename Scalar>
UnitAxis2D<Scalar>::UnitAxis2D(Scalar lambda, Scalar kappa)
{
  // Very small values are exactly zero
  lambda = std::abs(lambda) < eps_ ? 0 : lambda;
  kappa = std::abs(kappa) < eps_ ? 0 : kappa;
  // Set lambda = 1.0 if both are zero
  lambda = (std::abs(lambda) < eps_ && std::abs(kappa) < eps_) ? 1 : lambda;
  // Change signs if lambda is negative (order here is important!)
  kappa = lambda < 0 ? -kappa : kappa;
  lambda = lambda < 0 ? -lambda : lambda;
  // Ensure if lambda = 0, kappa = 1 (note lambda would be exactly 0 at this point)
  kappa = lambda == 0 ? 1 : kappa;
  // Normalize parameters (if required)
  Scalar square_sum = lambda * lambda + kappa * kappa;
  if (std::abs(square_sum - 1) > eps_)
  {
    Scalar normalizer = std::sqrt(square_sum);
    lambda_ = lambda / normalizer;
    kappa_ = kappa / normalizer;
  }
  else
  {
    lambda_ = lambda;
    kappa_ = kappa;
  }
}

// Operators

template <typename Scalar>
Scalar UnitAxis2D<Scalar>::operator[](std::size_t i) const
{
  if (i != 0 && i != 1)
  {
    throw std::out_of_range("UnitAxis2D<Scalar>::operator[]");
  }
  return i == 0 ? lambda_ : kappa_;
}

template <typename Scalar>
typename UnitAxis2D<Scalar>::Matrix UnitAxis2D<Scalar>::operator+() const
{
  Matrix ret;
  ret << lambda_, -kappa_, kappa_, lambda_;
  return ret;
}

template <typename Scalar>
typename UnitAxis2D<Scalar>::Matrix UnitAxis2D<Scalar>::operator-() const
{
  Matrix ret;
  ret << lambda_, kappa_, -kappa_, lambda_;
  return ret;
}

template <typename Scalar>
UnitAxis2D<Scalar>& UnitAxis2D<Scalar>::operator=(const UnitAxis2D<Scalar>& rhs)
{
  lambda_ = rhs.lambda();
  kappa_ = rhs.kappa();
  return *this;
}

// Non-member operators 

template <typename Scalar>
std::ostream& operator<<(std::ostream& os, const UnitAxis2D<Scalar>& m)
{
  os << "lambda: " << m.lambda() << ", kappa: " << m.kappa();
  return os;
}

template <typename Scalar>
bool operator==(const UnitAxis2D<Scalar>& m, const UnitAxis2D<Scalar>& n)
{
  Scalar eps = std::numeric_limits<Scalar>::epsilon();
  return std::abs(m.lambda() - n.lambda()) < eps && std::abs(m.kappa() - n.kappa()) < eps;
}

template <typename Scalar>
bool operator!=(const UnitAxis2D<Scalar>& m, const UnitAxis2D<Scalar>& n)
{
  return !(m == n);
}

template <typename Scalar>
UnitAxis2D<Scalar> operator*(const typename UnitAxis2D<Scalar>::Matrix& mat, const UnitAxis2D<Scalar>& n)
{
  return UnitAxis2D<Scalar>(mat * n.vector());
}

template <typename Scalar>
Scalar boxminus(const UnitAxis2D<Scalar>& m, const UnitAxis2D<Scalar>& n)
{
  return (-n * m).log();
}

// Non-member functions

template <typename Scalar>
UnitAxis2D<Scalar> boxplus(const UnitAxis2D<Scalar>& m, Scalar phi)
{
  return UnitAxis2D<Scalar>(+m * UnitAxis2D<Scalar>(phi));
}

template <typename Scalar>
Scalar distance(const UnitAxis2D<Scalar>& m, const UnitAxis2D<Scalar>& n)
{
  return std::acos(dot(m, n));
}

template <typename Scalar>
Scalar dot(const UnitAxis2D<Scalar>& m, const UnitAxis2D<Scalar>& n)
{
  Scalar x = std::abs(m.lambda() * n.lambda() + m.kappa() * n.kappa());
  return x > 1 ? 1 : x; // numerical check to ensure x <= 1.0
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis3D
///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
class UnitAxis3D
{
public:
  // Typedefs
  typedef Eigen::Matrix<Scalar, 2, 1> Vector2;
  typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
  typedef Eigen::Matrix<Scalar, 3, 3> Matrix;
  // Constructors
  UnitAxis3D(Scalar, Scalar, Scalar);
  UnitAxis3D(const Vector2&);
  UnitAxis3D() : UnitAxis3D(1, 0, 0) {}
  UnitAxis3D(const Vector3& v) : UnitAxis3D(v(0), v(1), v(2)) {}
  UnitAxis3D(Scalar lambda, const Vector2& kappa) : UnitAxis3D(lambda, kappa(0), kappa(1)) {}
  UnitAxis3D(const UnitAxis3D& orig) : lambda_(orig.lambda_), kappa_(orig.kappa_) {}
  // Instrospection
  Vector2     direction() const {return kappa_.norm() > eps_ ? kappa_.normalized() : Vector2(1, 0);}
  Scalar      inclination() const {return std::acos(lambda_);}
  Vector2     kappa() const {return kappa_;}
  Scalar      lambda() const {return lambda_;}
  Vector3     vector() const {return Vector3(lambda_, kappa_(0), kappa_(1));}
   // Mathematical methods
  UnitAxis3D  inv() const {return UnitAxis3D<Scalar>(lambda_, -kappa_);}
  Vector2     log() const {return direction() * inclination();}
  // Operators
  Scalar      operator[](std::size_t) const;
  Matrix      operator+() const;
  Matrix      operator-() const;
  UnitAxis3D& operator=(const UnitAxis3D&);
private:
  Scalar      lambda_;
  Vector2     kappa_;
  Scalar      eps_ = std::numeric_limits<Scalar>::epsilon();
  void        put_in_H2_(Scalar&, Vector2&);
};

// Constructors

template <typename Scalar>
UnitAxis3D<Scalar>::UnitAxis3D(Scalar lambda, Scalar kappa_1, Scalar kappa_2)
{
  Vector2 kappa(kappa_1, kappa_2);
  put_in_H2_(lambda, kappa);
  lambda_ = lambda;
  kappa_ = kappa;
}

template <typename Scalar>
UnitAxis3D<Scalar>::UnitAxis3D(const typename UnitAxis3D<Scalar>::Vector2& phi)
{
  Scalar lambda = 1;
  Vector2 kappa(0, 0);
  Scalar mag = phi.norm();
  if (mag > std::numeric_limits<Scalar>::epsilon())
  {
    lambda = std::cos(mag);
    kappa = std::sin(mag) * phi / mag;
  }
  put_in_H2_(lambda, kappa);
  lambda_ = lambda;
  kappa_ = kappa;
}

// Operators

template <typename Scalar>
Scalar UnitAxis3D<Scalar>::operator[](std::size_t i) const
{
  if (i != 0 && i != 1 && i != 2)
  {
    throw std::out_of_range("UnitAxis3D<Scalar>::operator[]");
  }
  return i == 0 ? lambda_ : kappa_(i - 1);
}

template <typename Scalar>
typename UnitAxis3D<Scalar>::Matrix UnitAxis3D<Scalar>::operator+() const
{
  Matrix ret;
  ret << lambda_, -kappa_(0), -kappa_(1), 
         kappa_(0), lambda_ + 1/(lambda_+1) * kappa_(1) * kappa_(1), -1/(lambda_+1) * kappa_(0) * kappa_(1),
         kappa_(1), -1/(lambda_+1) * kappa_(0) * kappa_(1), lambda_ + 1/(lambda_+1) * kappa_(0) * kappa_(0);
  return ret;
}

template <typename Scalar>
typename UnitAxis3D<Scalar>::Matrix UnitAxis3D<Scalar>::operator-() const
{
  return +(this->inv());
}

template <typename Scalar>
UnitAxis3D<Scalar>& UnitAxis3D<Scalar>::operator=(const UnitAxis3D<Scalar>& rhs)
{
  lambda_ = rhs.lambda();
  kappa_ = rhs.kappa();
  return *this;
}

// Utility methods

template <typename Scalar>
void UnitAxis3D<Scalar>::put_in_H2_(Scalar& lambda, Vector2& kappa)
{
  // Set to exactly zero if close to zero
  lambda = std::abs(lambda) < eps_ ? 0 : lambda;
  if (kappa.norm() < 0)
  {
    kappa.setZero();
  }

  // Set lambda = 1.0 if everything is zero
  lambda = (std::abs(lambda) < eps_ && kappa.norm() < eps_) ? 1 : lambda;

  // Change signs if lambda is negative
  if (lambda < 0)
  {
    kappa = -kappa;
    lambda = -lambda;
  }

  // Ensure if lambda = 0 and kappa_1 = 0, then kappa_2 = 1
  kappa(1) = (lambda == 0 && kappa(0) == 0) ? 1 : kappa(1);

  // Normalize parameters (if required)
  Scalar square_sum = lambda * lambda + kappa.squaredNorm();
  if (std::abs(square_sum - 1) > eps_)
  {
    Scalar normalizer = std::sqrt(square_sum);
    lambda = lambda / normalizer;
    kappa = kappa / normalizer;
  }
}

// Non-member operators

template <typename Scalar>
std::ostream& operator<<(std::ostream& os, const UnitAxis3D<Scalar>& m)
{
  os << "lambda: " << m.lambda() << ", kappa: [" << m.kappa()[0] << ", " << m.kappa()[1] << "]";
  return os;
}

template <typename Scalar>
bool operator==(const UnitAxis3D<Scalar>& m, const UnitAxis3D<Scalar>& n)
{
  Scalar eps = std::numeric_limits<Scalar>::epsilon();
  return std::abs(m.lambda() - n.lambda()) < eps && 
         std::abs(m.kappa()[0] - n.kappa()[0]) < eps && 
         std::abs(m.kappa()[1] - n.kappa()[1]) < eps;
}

template <typename Scalar>
bool operator!=(const UnitAxis3D<Scalar>& m, const UnitAxis3D<Scalar>& n)
{
  return !(m == n);
}

template <typename Scalar>
UnitAxis3D<Scalar> operator*(const typename UnitAxis3D<Scalar>::Matrix& mat, const UnitAxis3D<Scalar>& n)
{
  return UnitAxis3D<Scalar>(mat * n.vector());
}

template <typename Scalar>
typename UnitAxis3D<Scalar>::Vector2 boxminus(const UnitAxis3D<Scalar>& m, const UnitAxis3D<Scalar>& n)
{
  return (-n * m).log();
}

// Non-member functions

template <typename Scalar>
UnitAxis3D<Scalar> boxplus(const UnitAxis3D<Scalar>& m, const typename UnitAxis3D<Scalar>::Vector2& phi)
{
  return UnitAxis3D<Scalar>(+m * UnitAxis3D<Scalar>(phi));
}

template <typename Scalar>
Scalar distance(const UnitAxis3D<Scalar>& m, const UnitAxis3D<Scalar>& n)
{
  return std::acos(dot(m, n));
}

template <typename Scalar>
Scalar dot(const UnitAxis3D<Scalar>& m, const UnitAxis3D<Scalar>& n)
{
  Scalar x = std::abs(m.lambda() * n.lambda() + m.kappa().dot(n.kappa));
  return x > 1 ? 1 : x; // numerical check to ensure x <= 1.0
}

} // namespace axes

#endif
