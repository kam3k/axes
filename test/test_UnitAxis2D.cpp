#include "catch.h"
#include "axes.h"
#include <iostream>

using namespace axes;

TEST_CASE("2D constructors")
{

  SECTION("Default constructor")
  {
    UnitAxis2D<double> m;
    REQUIRE(m.lambda() == 1.0);
    REQUIRE(m.kappa() == 0.0);
    UnitAxis2D<long> n;
    REQUIRE(n.lambda() == 1.0L);
    REQUIRE(n.kappa() == 0.0L);
  }

  SECTION("Two parameter constructor")
  {
    UnitAxis2D<double> a(1.0, 0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa() == 0.0);

    UnitAxis2D<double> b(2.0, 0.0);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa() == 0.0);

    UnitAxis2D<double> c(1.0, 1.0);
    REQUIRE(c.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(c.kappa() == Approx(std::sqrt(2.0) / 2.0));

    UnitAxis2D<double> d(-1.0, 1.0);
    REQUIRE(d.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(d.kappa() == Approx(-std::sqrt(2.0) / 2.0));

    UnitAxis2D<double> e(1.0, -1.0);
    REQUIRE(e.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(e.kappa() == Approx(-std::sqrt(2.0) / 2.0));

    UnitAxis2D<double> f(0.0, 0.0);
    REQUIRE(f.lambda() == 1.0);
    REQUIRE(f.kappa() == 0.0);

    UnitAxis2D<double> g(0.0, -1.0);
    REQUIRE(g.lambda() == 0.0);
    REQUIRE(g.kappa() == 1.0);
  }

  SECTION("Angle constructor")
  {
    const long double pi = 3.141592653589793238462643383279502884L;

    UnitAxis2D<double> a(0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa() == 0.0);

    UnitAxis2D<double> b(pi);
    REQUIRE(b.lambda() == Approx(1.0));
    REQUIRE(b.kappa() == Approx(0.0));

    UnitAxis2D<double> c(pi/2.0);
    REQUIRE(c.lambda() == Approx(0.0));
    REQUIRE(c.kappa() == Approx(1.0));

    UnitAxis2D<double> d(pi/3.0);
    REQUIRE(d.lambda() == Approx(0.5));
    REQUIRE(d.kappa() == Approx(std::sqrt(3.0) / 2.0));

    UnitAxis2D<double> e(3.0*pi/4.0);
    REQUIRE(e.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(e.kappa() == Approx(-std::sqrt(2.0) / 2.0));
  }

  SECTION("2-vector constructor")
  {
    Eigen::Vector2d u(1.0, 1.0);
    UnitAxis2D<double> a(u);
    REQUIRE(a.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(a.kappa() == Approx(std::sqrt(2.0) / 2.0));

    Eigen::Vector2d v(0.0, 0.0);
    UnitAxis2D<double> b(v);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa() == 0.0);
  }

  SECTION("Copy constructor")
  {
    UnitAxis2D<double> a(1.0, 1.0);
    UnitAxis2D<double> b(a);
    REQUIRE(a.lambda() == b.lambda());
    REQUIRE(a.kappa() == b.kappa());

    UnitAxis2D<double> c(-1.0, 4.32);
    UnitAxis2D<double> d(c);
    REQUIRE(c.lambda() == d.lambda());
    REQUIRE(c.kappa() == d.kappa());
  }

}

TEST_CASE("2D introspection")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  UnitAxis2D<double> m;
  UnitAxis2D<double> n(1.0, 1.0);
  UnitAxis2D<double> p(0.0, -1.0);

  SECTION("Angle")
  {
    REQUIRE(m.angle() == Approx(0.0));
    REQUIRE(n.angle() == Approx(pi/4.0));
    REQUIRE(p.angle() == Approx(pi/2.0));
  }

  SECTION("Kappa")
  {
    REQUIRE(m.kappa() == 0.0);
    REQUIRE(n.kappa() == Approx(std::sqrt(2.0)/2));
    REQUIRE(p.kappa() == 1.0);
  }

  SECTION("Lambda")
  {
    REQUIRE(m.lambda() == 1.0);
    REQUIRE(n.lambda() == Approx(std::sqrt(2.0)/2));
    REQUIRE(p.lambda() == 0.0);
  }

  SECTION("Vector")
  {
    Eigen::Vector2d m_vec(1.0, 0.0);
    Eigen::Vector2d n_vec(std::sqrt(2.0)/2, std::sqrt(2.0)/2);
    Eigen::Vector2d p_vec(0.0, 1.0);

    REQUIRE(m.vector() == m_vec);
    REQUIRE(n.vector()[0] == Approx(n_vec[0]));
    REQUIRE(n.vector()[1] == Approx(n_vec[1]));
    REQUIRE(p.vector() == p_vec);
  }
}

TEST_CASE("2D mathematical methods")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  UnitAxis2D<double> m;
  UnitAxis2D<double> n(1.0, 1.0);
  UnitAxis2D<double> p(0.0, -1.0);
  UnitAxis2D<double> q(2.3, -3.7);

  SECTION("Inv")
  {
    REQUIRE(m.inv() == UnitAxis2D<double>(1.0, 0.0));
    REQUIRE(n.inv() == UnitAxis2D<double>(1.0, -1.0));
    REQUIRE(p.inv() == UnitAxis2D<double>(0.0, 1.0));
    REQUIRE(q.inv() == UnitAxis2D<double>(2.3, 3.7));
  }

  SECTION("Log")
  {
    REQUIRE(m.log() == Approx(0.0));
    REQUIRE(n.log() == Approx(pi/4.0));
    REQUIRE(p.log() == Approx(pi/2.0));
    REQUIRE(q.log() == Approx(-1.014630096674437));
  }
}

TEST_CASE("2D operators")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  UnitAxis2D<double> m;
  UnitAxis2D<double> n(1.0, 1.0);
  UnitAxis2D<double> p(0.0, -1.0);

  SECTION("Indexing")
  {
    REQUIRE(m[0] == 1.0);
    REQUIRE(m[1] == 0.0);
    REQUIRE(n[0] == Approx(std::sqrt(2.0)/2));
    REQUIRE(n[1] == Approx(std::sqrt(2.0)/2));
    REQUIRE(p[0] == 0.0);
    REQUIRE(p[1] == 1.0);
  }

  SECTION("+ compound")
  {
    Eigen::Matrix2d m_plus;
    m_plus << m.lambda(), -m.kappa(), m.kappa(), m.lambda();
    Eigen::Matrix2d n_plus;
    n_plus << n.lambda(), -n.kappa(), n.kappa(), n.lambda();
    Eigen::Matrix2d p_plus;
    p_plus << p.lambda(), -p.kappa(), p.kappa(), p.lambda();

    Eigen::Matrix2d m_plus_calc = +m;
    Eigen::Matrix2d n_plus_calc = +n;
    Eigen::Matrix2d p_plus_calc = +p;
    REQUIRE(m_plus_calc(0, 0) == Approx(m_plus(0, 0)));
    REQUIRE(m_plus_calc(1, 0) == Approx(m_plus(1, 0)));
    REQUIRE(m_plus_calc(0, 1) == Approx(m_plus(0, 1)));
    REQUIRE(m_plus_calc(1, 1) == Approx(m_plus(1, 1)));
    REQUIRE(n_plus_calc(0, 0) == Approx(n_plus(0, 0)));
    REQUIRE(n_plus_calc(1, 0) == Approx(n_plus(1, 0)));
    REQUIRE(n_plus_calc(0, 1) == Approx(n_plus(0, 1)));
    REQUIRE(n_plus_calc(1, 1) == Approx(n_plus(1, 1)));
    REQUIRE(p_plus_calc(0, 0) == Approx(p_plus(0, 0)));
    REQUIRE(p_plus_calc(1, 0) == Approx(p_plus(1, 0)));
    REQUIRE(p_plus_calc(0, 1) == Approx(p_plus(0, 1)));
    REQUIRE(p_plus_calc(1, 1) == Approx(p_plus(1, 1)));
  }

  SECTION("- compound")
  {
    Eigen::Matrix2d m_minus;
    m_minus << m.lambda(), m.kappa(), -m.kappa(), m.lambda();
    Eigen::Matrix2d n_minus;
    n_minus << n.lambda(), n.kappa(), -n.kappa(), n.lambda();
    Eigen::Matrix2d p_minus;
    p_minus << p.lambda(), p.kappa(), -p.kappa(), p.lambda();

    Eigen::Matrix2d m_minus_calc = -m;
    Eigen::Matrix2d n_minus_calc = -n;
    Eigen::Matrix2d p_minus_calc = -p;
    REQUIRE(m_minus_calc(0, 0) == Approx(m_minus(0, 0)));
    REQUIRE(m_minus_calc(1, 0) == Approx(m_minus(1, 0)));
    REQUIRE(m_minus_calc(0, 1) == Approx(m_minus(0, 1)));
    REQUIRE(m_minus_calc(1, 1) == Approx(m_minus(1, 1)));
    REQUIRE(n_minus_calc(0, 0) == Approx(n_minus(0, 0)));
    REQUIRE(n_minus_calc(1, 0) == Approx(n_minus(1, 0)));
    REQUIRE(n_minus_calc(0, 1) == Approx(n_minus(0, 1)));
    REQUIRE(n_minus_calc(1, 1) == Approx(n_minus(1, 1)));
    REQUIRE(p_minus_calc(0, 0) == Approx(p_minus(0, 0)));
    REQUIRE(p_minus_calc(1, 0) == Approx(p_minus(1, 0)));
    REQUIRE(p_minus_calc(0, 1) == Approx(p_minus(0, 1)));
    REQUIRE(p_minus_calc(1, 1) == Approx(p_minus(1, 1)));
  }

  SECTION("Assignment")
  {
    m = UnitAxis2D<double>(4.3, -1.1);
    n = UnitAxis2D<double>(-1.3, -11.9);
    p = UnitAxis2D<double>(7.2, 4.8);

    REQUIRE(m == UnitAxis2D<double>(4.3, -1.1));
    REQUIRE(n == UnitAxis2D<double>(-1.3, -11.9));
    REQUIRE(p == UnitAxis2D<double>(7.2, 4.8));
  }
}

TEST_CASE("2D functions")
{
  const double pi = 3.141592653589793238462643383279502884;
  UnitAxis2D<double> m;
  UnitAxis2D<double> n(1.0, 1.0);
  UnitAxis2D<double> p(0.0, -1.0);
  UnitAxis2D<double> q(std::sqrt(3.0)/2, 0.5);
  UnitAxis2D<double> r(0.5, -std::sqrt(3.0)/2);

  SECTION("Equivalency operators")
  {
    UnitAxis2D<double> s;
    UnitAxis2D<double> t(1.0, 1.0);
    UnitAxis2D<double> u(0.0, -1.0);

    REQUIRE(m == m);
    REQUIRE(m == s);
    REQUIRE(n == n);
    REQUIRE(n == t);
    REQUIRE(p == p);
    REQUIRE(p == u);

    REQUIRE(m != n);
    REQUIRE(m != p);
    REQUIRE(m != t);
    REQUIRE(m != u);
  }

  SECTION("Axis product")
  {
    REQUIRE(+m * UnitAxis2D<double>() == m); 
    REQUIRE(+n * UnitAxis2D<double>() == n); 
    REQUIRE(+p * UnitAxis2D<double>() == p); 
    REQUIRE(+q * q == UnitAxis2D<double>(0.5, std::sqrt(3.0)/2));
    REQUIRE(+q * r == UnitAxis2D<double>(std::sqrt(3.0)/2, -0.5));
    REQUIRE(+r * q == UnitAxis2D<double>(std::sqrt(3.0)/2, -0.5));
    REQUIRE(+q * +r * q == UnitAxis2D<double>());
    REQUIRE(+r * +q * q == UnitAxis2D<double>());
    REQUIRE(+r * r == UnitAxis2D<double>(0.5, std::sqrt(3.0)/2));

    REQUIRE(-m * UnitAxis2D<double>() == m.inv()); 
    REQUIRE(-n * UnitAxis2D<double>() == n.inv()); 
    REQUIRE(-p * UnitAxis2D<double>() == p.inv()); 
    REQUIRE(-q * q == UnitAxis2D<double>());
    REQUIRE(-q * r == UnitAxis2D<double>(0, 1));
    REQUIRE(-r * q == UnitAxis2D<double>(0, 1));
    REQUIRE(-q * -r * q == UnitAxis2D<double>(0.5, std::sqrt(3.0)/2));
    REQUIRE(-r * -q * q == UnitAxis2D<double>(0.5, std::sqrt(3.0)/2));
    REQUIRE(-r * r == UnitAxis2D<double>());
  }

  SECTION("Boxminus")
  {
    REQUIRE(boxminus(m, m) == Approx(0));
    REQUIRE(boxminus(n, n) == Approx(0));
    REQUIRE(boxminus(q, r) == Approx(pi/2));
    REQUIRE(boxminus(r, q) == Approx(pi/2));
    REQUIRE(boxminus(p, q) == Approx(pi/3));
    REQUIRE(boxminus(q, p) == Approx(-pi/3));
    REQUIRE(boxminus(n, r) == Approx(7*pi/12 - pi));
    REQUIRE(boxminus(r, n) == Approx(-7*pi/12 + pi));
  }

  SECTION("Boxplus")
  {
    REQUIRE(boxplus(m, 0.0) == m);
    REQUIRE(boxplus(m, pi/4) == n);
    REQUIRE(boxplus(m, -pi/3) == r);
    REQUIRE(boxplus(q, -pi/2) == r);
    REQUIRE(boxplus(r, pi/2) == q);
  }

  SECTION("Manifold axioms")
  {
    REQUIRE(boxplus(m, 0.0) == m);
    REQUIRE(boxplus(r, boxminus(q, r)) == q);
    REQUIRE(boxplus(n, boxminus(r, n)) == r);
    REQUIRE(boxminus(boxplus(q, 0.55), q) == Approx(0.55));
    REQUIRE(boxminus(boxplus(r, -0.32), r) == Approx(-0.32));
  }

  SECTION("Distance")
  {
    REQUIRE(distance(m, n) == Approx(pi/4));
    REQUIRE(distance(n, m) == Approx(pi/4));
    REQUIRE(distance(q, r) == Approx(pi/2));
    REQUIRE(distance(r, q) == Approx(pi/2));
    REQUIRE(distance(m, m) == Approx(0));
    REQUIRE(distance(n, n) == Approx(0));
    REQUIRE(distance(p, q) == Approx(pi/3));
    REQUIRE(distance(q, p) == Approx(pi/3));
    REQUIRE(distance(n, r) == Approx(pi - 7*pi/12));
    REQUIRE(distance(r, n) == Approx(pi - 7*pi/12));
  }

  SECTION("Dot")
  {
    REQUIRE(dot(m, n) == Approx(std::cos(pi/4)));
    REQUIRE(dot(n, m) == Approx(std::cos(pi/4)));
    REQUIRE(dot(q, r) == Approx(std::cos(pi/2)));
    REQUIRE(dot(r, q) == Approx(std::cos(pi/2)));
    REQUIRE(dot(m, m) == Approx(std::cos(0)));
    REQUIRE(dot(n, n) == Approx(std::cos(0)));
    REQUIRE(dot(p, q) == Approx(std::cos(pi/3)));
    REQUIRE(dot(q, p) == Approx(std::cos(pi/3)));
    REQUIRE(dot(n, r) == Approx(std::cos(pi - 7*pi/12)));
    REQUIRE(dot(r, n) == Approx(std::cos(pi - 7*pi/12)));
  }
}
