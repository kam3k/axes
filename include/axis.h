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

template <typename T>
class UnitAxis2D
{
public:
  // Typedefs
  typedef Eigen::Matrix<T, 2, 1> Vector;
  typedef Eigen::Matrix<T, 2, 2> Matrix;
  // Constructors
  UnitAxis2D(T, T);
  UnitAxis2D(T phi) : UnitAxis2D(std::cos(phi), std::sin(phi)) {}
  UnitAxis2D() : UnitAxis2D(1, 0) {}
  UnitAxis2D(const Vector v) : UnitAxis2D(v(0), v(1)) {}
  UnitAxis2D(const UnitAxis2D& orig) : lambda_(orig.lambda_), kappa_(orig.kappa_) {}
  // Introspection
  T           angle() const {return std::atan(kappa_ / lambda_);}
  T           kappa() const {return kappa_;}
  T           lambda() const {return lambda_;}
  Vector      vector() const {return Vector(lambda_, kappa_);}
  // Mathematical methods
  UnitAxis2D  inv() const {return UnitAxis2D<T>(lambda_, -kappa_);}
  T           log() const {return angle();}
  // Operators
  T           operator[](std::size_t) const;
  Matrix      operator+() const;
  Matrix      operator-() const;
  UnitAxis2D& operator=(const UnitAxis2D&);
private:
  T           lambda_;
  T           kappa_;
  T           eps_ = std::numeric_limits<T>::epsilon();
};

// Constructor

template <typename T>
UnitAxis2D<T>::UnitAxis2D(T lambda, T kappa)
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
  double square_sum = lambda * lambda + kappa * kappa;
  if (std::abs(square_sum - 1) > eps_)
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

// Operators

template <typename T>
T UnitAxis2D<T>::operator[](std::size_t i) const
{
  if (i != 0 && i != 1)
  {
    throw std::out_of_range("UnitAxis2D<T>::operator[]");
  }
  return i == 0 ? lambda_ : kappa_;
}

template <typename T>
typename UnitAxis2D<T>::Matrix UnitAxis2D<T>::operator+() const
{
  Matrix ret;
  ret << lambda_, -kappa_, kappa_, lambda_;
  return ret;
}

template <typename T>
typename UnitAxis2D<T>::Matrix UnitAxis2D<T>::operator-() const
{
  Matrix ret;
  ret << lambda_, kappa_, -kappa_, lambda_;
  return ret;
}

template <typename T>
UnitAxis2D<T>& UnitAxis2D<T>::operator=(const UnitAxis2D<T>& rhs)
{
  lambda_ = rhs.lambda();
  kappa_ = rhs.kappa();
  return *this;
}

// Non-member operators 

template <typename T>
std::ostream& operator<<(std::ostream& os, const UnitAxis2D<T>& m)
{
  os << "lambda: " << m.lambda() << ", kappa: " << m.kappa();
  return os;
}

template <typename T>
bool operator==(const UnitAxis2D<T>& m, const UnitAxis2D<T>& n)
{
  T eps = std::numeric_limits<T>::epsilon();
  return std::abs(m.lambda() - n.lambda()) < eps && std::abs(m.kappa() - n.kappa()) < eps;
}

template <typename T>
bool operator!=(const UnitAxis2D<T>& m, const UnitAxis2D<T>& n)
{
  return !(m == n);
}

template <typename T>
UnitAxis2D<T> operator*(const typename UnitAxis2D<T>::Matrix& mat, const UnitAxis2D<T>& n)
{
  return UnitAxis2D<T>(mat * n.vector());
}

template <typename T>
double boxminus(const UnitAxis2D<T>& m, const UnitAxis2D<T>& n)
{
  return (-n * m).log();
}

// Non-member functions

template <typename T>
UnitAxis2D<T> boxplus(const UnitAxis2D<T>& m, double phi)
{
  return UnitAxis2D<T>(+m * UnitAxis2D<T>(phi));
}

template <typename T>
T distance(const UnitAxis2D<T>& m, const UnitAxis2D<T>& n)
{
  return std::acos(dot(m, n));
}

template <typename T>
T dot(const UnitAxis2D<T>& m, const UnitAxis2D<T>& n)
{
  double x = std::abs(m.lambda() * n.lambda() + m.kappa() * n.kappa());
  return x > 1 ? 1 : x; // numerical check to ensure x <= 1.0
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// UnitAxis3D
///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
class UnitAxis3D
{
public:
  // Typedegs
  typedef Eigen::Matrix<T, 2, 1> Vector2;
  typedef Eigen::Matrix<T, 3, 1> Vector3;
  typedef Eigen::Matrix<T, 3, 3> Matrix;
  // Constructors
  UnitAxis3D(T, T, T);
  UnitAxis3D(const Vector2&);
  UnitAxis3D() : UnitAxis3D(1, 0, 0) {}
  UnitAxis3D(const Vector3& v) : UnitAxis3D(v(0), v(1), v(2)) {}
  UnitAxis3D(double lambda, const Vector2& kappa) : UnitAxis3D(lambda, kappa(0), kappa(1)) {}
  UnitAxis3D(const UnitAxis3D& orig) : lambda_(orig.lambda_), kappa_(orig.kappa_) {}
  // Instrospection
  Vector2     direction() const {return kappa_.normalized();}
  T           inclination() const {return std::acos(lambda_);}
  Vector2     kappa() const {return kappa_;}
  T           lambda() const {return lambda_;}
  Vector3     vector() const {return Vector3(lambda_, kappa_(0), kappa_(1));}
   // Mathematical methods
  UnitAxis3D  inv() const {return UnitAxis3D<T>(lambda_, -kappa_);}
  Vector2     log() const {return direction() * inclination();}
  // Operators
  T           operator[](std::size_t) const;
  Matrix      operator+() const;
  Matrix      operator-() const;
  UnitAxis3D& operator=(const UnitAxis3D&);
private:
  T           lambda_;
  Vector2     kappa_;
  T           eps_ = std::numeric_limits<T>::epsilon();
};

// Constructors

template <typename T>
UnitAxis3D<T>::UnitAxis3D(T lambda, T kappa_1, T kappa_2)
{
  typename UnitAxis3D<T>::Vector2 kappa(kappa_1, kappa_2);
  // Very small values are exactly zero
  lambda = std::abs(lambda) < eps_ ? 0 : lambda;
  if (kappa.norm() < 0)
  {
    kappa.setZero();
  }
  // Set lambda = 1.0 if everything is zero
  lambda = (std::abs(lambda) < eps_ && kappa.norm() < eps_) ? 1 : lambda;
  // Change signs if lambda is negative (order here is important!)
  kappa = lambda < 0 ? -kappa : kappa;
  lambda = lambda < 0 ? -lambda : lambda;
  // Ensure if lambda = 0 and kappa_1 = 0, then kappa_2 = 1
  kappa(1) = (lambda == 0 && kappa(0) == 0) ? 1 : kappa(1);
  // Normalize parameters (if required)
  double square_sum = lambda * lambda + kappa.squaredNorm();
  if (std::abs(square_sum - 1) > eps_)
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

template <typename T>
UnitAxis3D<T>::UnitAxis3D(const typename UnitAxis3D<T>::Vector2& phi)
{
  double mag = phi.norm();
  if (mag < std::numeric_limits<T>::epsilon())
  {
    return UnitAxis3D<T>(1, 0, 0);
  }
  typename UnitAxis3D<T>::Vector2 kappa = std::sin(mag) * phi / mag;
  return UnitAxis3D<T>(std::cos(mag), kappa(0), kappa(1));
}

// Operators

template <typename T>
T UnitAxis3D<T>::operator[](std::size_t i) const
{
  if (i != 0 && i != 1 && i != 2)
  {
    throw std::out_of_range("UnitAxis3D<T>::operator[]");
  }
  return i == 0 ? lambda_ : kappa_(i - 1);
}

template <typename T>
typename UnitAxis3D<T>::Matrix UnitAxis3D<T>::operator+() const
{
  Matrix ret;
  ret << lambda_, -kappa_(0), -kappa_(1), 
         kappa(0), lambda_ + 1/(lambda_+1) * kappa(1) * kappa(1), -1/(kappa(0) * kappa(1)),
         kappa(1), lambda_ + 1/(lambda_+1) * kappa(0) * kappa(0), -1/(kappa(0) * kappa(1));
  return ret;
}

template <typename T>
typename UnitAxis3D<T>::Matrix UnitAxis3D<T>::operator-() const
{
  return +(this->inv());
}

template <typename T>
UnitAxis3D<T>& UnitAxis3D<T>::operator=(const UnitAxis3D<T>& rhs)
{
  lambda_ = rhs.lambda();
  kappa_ = rhs.kappa();
  return *this;
}

// Non-member operators

// Functions acting on UnitAxis3D
//std::ostream    &operator<<(std::ostream&, const UnitAxis3D&);
//bool            operator==(const UnitAxis3D&, const UnitAxis3D&);
//bool            operator!=(const UnitAxis3D&, const UnitAxis3D&);
//UnitAxis3D      operator+(const UnitAxis3D&, const UnitAxis3D&);
//Eigen::Vector2d boxminus(const UnitAxis3D&, const UnitAxis3D&);
//UnitAxis3D      boxplus(const UnitAxis3D&, const Eigen::Vector3d&);
//double          distance(const UnitAxis3D&, const UnitAxis3D&);
//double          dot(const UnitAxis3D&, const UnitAxis3D&);

} // namespace axes

#endif
