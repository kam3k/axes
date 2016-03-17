#include "catch.h"
#include "axes.h"
#include <iostream>

using namespace axes;

TEST_CASE("3D constructors")
{

  SECTION("Default constructor")
  {
    UnitAxis3D<double> m;
    REQUIRE(m.lambda() == 1.0);
    REQUIRE(m.kappa()[0] == 0.0);
    REQUIRE(m.kappa()[1] == 0.0);
  }

  SECTION("Three parameter constructor")
  {
    UnitAxis3D<double> a(1.0, 0.0, 0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa()[0] == 0.0);
    REQUIRE(a.kappa()[1] == 0.0);

    UnitAxis3D<double> b(2.0, 0.0, 0.0);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa()[0] == 0.0);
    REQUIRE(b.kappa()[1] == 0.0);

    UnitAxis3D<double> c(1.0, 1.0, 0.0);
    REQUIRE(c.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(c.kappa()[0] == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(c.kappa()[1] == 0.0);

    UnitAxis3D<double> d(-1.0, 1.0, 0.0);
    REQUIRE(d.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(d.kappa()[0] == Approx(-std::sqrt(2.0) / 2.0));
    REQUIRE(d.kappa()[1] == 0.0);

    UnitAxis3D<double> e(0.0, 0.0, -1.0);
    REQUIRE(e.lambda() == 0.0);
    REQUIRE(e.kappa()[0] == 0.0);
    REQUIRE(e.kappa()[1] == 1.0);

    UnitAxis3D<double> f(0.0, 0.0, 0.0);
    REQUIRE(f.lambda() == 1.0);
    REQUIRE(f.kappa()[0] == 0.0);
    REQUIRE(f.kappa()[1] == 0.0);

    UnitAxis3D<double> g(-1.0, -1.0, -1.0);
    REQUIRE(g.lambda() == Approx(1 / std::sqrt(3.0)));
    REQUIRE(g.kappa()[0] == Approx(1 / std::sqrt(3.0)));
    REQUIRE(g.kappa()[1] == Approx(1 / std::sqrt(3.0)));
  }

  SECTION("Axis vector constructor")
  {
    const long double pi = 3.141592653589793238462643383279502884L;
    Eigen::Vector2d e1(1, 0);
    Eigen::Vector2d e2(0, 1);

    UnitAxis3D<double> a(Eigen::Vector2d(0, 0));
    REQUIRE(a == UnitAxis3D<double>());

    UnitAxis3D<double> b(Eigen::Vector2d(pi/2, 0));
    REQUIRE(b.lambda() == 0.0);
    REQUIRE(b.kappa()[0] == 1.0);
    REQUIRE(b.kappa()[1] == 0.0);

    UnitAxis3D<double> c(Eigen::Vector2d(0, -pi/2));
    REQUIRE(c.lambda() == Approx(0.0));
    REQUIRE(c.kappa()[0] == Approx(0.0));
    REQUIRE(c.kappa()[1] == Approx(1.0));

    UnitAxis3D<double> d(Eigen::Vector2d(pi/3, -pi/3));
    REQUIRE(d.lambda() == Approx(0.08971456178219864));
    REQUIRE(d.kappa()[0] == Approx(0.704255));
    REQUIRE(d.kappa()[1] == Approx(-0.704255));

    UnitAxis3D<double> e(Eigen::Vector2d(-1.3, 2.4));
    REQUIRE(e.lambda() == Approx(0.9162721722720432));
    REQUIRE(e.kappa()[0] == Approx(0.19077819007809607));
    REQUIRE(e.kappa()[1] == Approx(-0.35220588937494657));
  }

  SECTION("3-vector constructor")
  {
    UnitAxis3D<double> a(Eigen::Vector3d(1.0, 1.0, 0.0));
    REQUIRE(a.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(a.kappa()[0] == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(a.kappa()[1] == 0.0);

    UnitAxis3D<double> b(Eigen::Vector3d(-1.0, 1.0, 0.0));
    REQUIRE(b.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(b.kappa()[0] == Approx(-std::sqrt(2.0) / 2.0));
    REQUIRE(b.kappa()[1] == 0.0);

    UnitAxis3D<double> c(Eigen::Vector3d(0.0, 0.0, -1.0));
    REQUIRE(c.lambda() == 0.0);
    REQUIRE(c.kappa()[0] == 0.0);
    REQUIRE(c.kappa()[1] == 1.0);

    UnitAxis3D<double> d(Eigen::Vector3d(0.0, 0.0, 0.0));
    REQUIRE(d.lambda() == 1.0);
    REQUIRE(d.kappa()[0] == 0.0);
    REQUIRE(d.kappa()[1] == 0.0);
  }

  SECTION("Copy constructor")
  {
    UnitAxis3D<double> a(1.0, 1.0, 1.0);
    UnitAxis3D<double> b(a);
    REQUIRE(a.lambda() == b.lambda());
    REQUIRE(a.kappa() == b.kappa());

    UnitAxis3D<double> c(-1.0, 4.32, -7.11);
    UnitAxis3D<double> d(c);
    REQUIRE(c.lambda() == d.lambda());
    REQUIRE(c.kappa() == d.kappa());
  }

}

TEST_CASE("3D introspection")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  UnitAxis3D<double> m;
  UnitAxis3D<double> n(1.0, 1.0, 0.0);
  UnitAxis3D<double> p(0.0, 0.0, -1.0);

  SECTION("Direction")
  {
    REQUIRE(m.direction() == Eigen::Vector2d(1, 0));
    REQUIRE(n.direction() == Eigen::Vector2d(1, 0));
    REQUIRE(p.direction() == Eigen::Vector2d(0, 1));
  }

  SECTION("Inclination")
  {
    REQUIRE(m.inclination() == Approx(0.0));
    REQUIRE(n.inclination() == Approx(pi/4.0));
    REQUIRE(p.inclination() == Approx(pi/2.0));
  }

  SECTION("Kappa")
  {
    REQUIRE(m.kappa() == Eigen::Vector2d(0, 0));
    REQUIRE(n.kappa()[0] == Approx(std::sqrt(2)/2));
    REQUIRE(n.kappa()[1] == 0.0);
    REQUIRE(p.kappa() == Eigen::Vector2d(0, 1));
  }

  SECTION("Lambda")
  {
    REQUIRE(m.lambda() == 1.0);
    REQUIRE(n.lambda() == Approx(std::sqrt(2.0)/2));
    REQUIRE(p.lambda() == 0.0);
  }

  SECTION("Vector")
  {
    Eigen::Vector3d m_vec(1.0, 0.0, 0.0);
    Eigen::Vector3d n_vec(std::sqrt(2)/2, std::sqrt(2)/2, 0.0);
    Eigen::Vector3d p_vec(0.0, 0.0, 1.0);

    REQUIRE(m.vector() == m_vec);
    REQUIRE(n.vector()[0] == Approx(n_vec[0]));
    REQUIRE(n.vector()[1] == Approx(n_vec[1]));
    REQUIRE(n.vector()[2] == Approx(n_vec[2]));
    REQUIRE(p.vector() == p_vec);
  }
}

TEST_CASE("3D mathematical methods")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  UnitAxis3D<double> m;
  UnitAxis3D<double> n(1.0, 1.0, 1.0);
  UnitAxis3D<double> p(0.0, 0.0, -1.0);
  UnitAxis3D<double> q(2.3, -3.7, -44.31);

  SECTION("Inv")
  {
    REQUIRE(m.inv() == UnitAxis3D<double>(1.0, 0.0, 0.0));
    REQUIRE(n.inv() == UnitAxis3D<double>(1.0, -1.0, -1.0));
    REQUIRE(p.inv() == UnitAxis3D<double>(0.0, 0.0, 1.0));
    REQUIRE(q.inv() == UnitAxis3D<double>(2.3, 3.7, 44.31));
  }

  SECTION("Log")
  {
    REQUIRE(m.log() == Eigen::Vector2d(0, 0));
    REQUIRE(n.log()[0] == Approx(0.67551085885604));
    REQUIRE(n.log()[1] == Approx(0.67551085885604));
    REQUIRE(p.log()[0] == 0.0);
    REQUIRE(p.log()[1] == Approx(pi/2));
    REQUIRE(q.log()[0] == Approx(-0.12641013466187973));
    REQUIRE(q.log()[1] == Approx(-1.5138467748291597));
  }
}

TEST_CASE("3D operators")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  UnitAxis3D<double> m;
  UnitAxis3D<double> n(1.0, 1.0, 0.0);
  UnitAxis3D<double> p(0.0, 0.0, -1.0);

  SECTION("Indexing")
  {
    REQUIRE(m[0] == 1.0);
    REQUIRE(m[1] == 0.0);
    REQUIRE(m[2] == 0.0);
    REQUIRE(n[0] == Approx(std::sqrt(2.0)/2));
    REQUIRE(n[1] == Approx(std::sqrt(2.0)/2));
    REQUIRE(n[2] == 0.0);
    REQUIRE(p[0] == 0.0);
    REQUIRE(p[1] == 0.0);
    REQUIRE(p[2] == 1.0);
  }

  SECTION("+ compound")
  {
    Eigen::Matrix3d m_plus;
    m_plus << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    Eigen::Matrix3d n_plus;
    n_plus << std::sqrt(2)/2, -std::sqrt(2)/2, 0, std::sqrt(2)/2, std::sqrt(2)/2, 0, 0, 0, 1;
    Eigen::Matrix3d p_plus;
    p_plus << 0, 0, 0, 0, 1, 0, 0, 0, 0;


    Eigen::Matrix3d m_plus_calc = +m;
    Eigen::Matrix3d n_plus_calc = +n;
    Eigen::Matrix3d p_plus_calc = +p;
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
    Eigen::Matrix3d m_minus;
    m_minus << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    Eigen::Matrix3d n_minus;
    n_minus << std::sqrt(2)/2, std::sqrt(2)/2, 0, -std::sqrt(2)/2, std::sqrt(2)/2, 0, 0, 0, 1;
    Eigen::Matrix3d p_minus;
    p_minus << 0, 0, 0, 0, 1, 0, 0, 0, 0;

    Eigen::Matrix3d m_minus_calc = -m;
    Eigen::Matrix3d n_minus_calc = -n;
    Eigen::Matrix3d p_minus_calc = -p;
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
    m = UnitAxis3D<double>(4.3, -1.1, 7.4);
    n = UnitAxis3D<double>(-1.3, -11.9, 1.1);
    p = UnitAxis3D<double>(7.2, 4.8, 9.9);

    REQUIRE(m == UnitAxis3D<double>(4.3, -1.1, 7.4));
    REQUIRE(n == UnitAxis3D<double>(-1.3, -11.9, 1.1));
    REQUIRE(p == UnitAxis3D<double>(7.2, 4.8, 9.9));
  }
}

TEST_CASE("3D functions")
{
  const double pi = 3.141592653589793238462643383279502884;
  UnitAxis3D<double> m;
  UnitAxis3D<double> n(1.0, 1.0, 0.0);
  UnitAxis3D<double> p(0.0, 0.0, -1.0);
  UnitAxis3D<double> q(std::sqrt(3.0)/2, 0.5, 0);
  UnitAxis3D<double> r(0.0, 0.5, -std::sqrt(3.0)/2);

  SECTION("Equivalency operators")
  {
    UnitAxis3D<double> s;
    UnitAxis3D<double> t(1.0, 1.0, 0.0);
    UnitAxis3D<double> u(0.0, 0.0, -1.0);

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
    REQUIRE(+m * UnitAxis3D<double>() == m); 
    REQUIRE(+n * UnitAxis3D<double>() == n); 
    REQUIRE(+p * UnitAxis3D<double>() == p); 
    REQUIRE(+q * q == UnitAxis3D<double>(0.5, std::sqrt(3.0)/2, 0));
    REQUIRE(+q * r == UnitAxis3D<double>(0.25, -0.43301270189221935, 0.8660254037844387));
    REQUIRE(+r * q == UnitAxis3D<double>(0.25, -0.8080127018922195, 0.5334936490538903));
    REQUIRE(+q * +r * q == UnitAxis3D<double>(0.6205127018922194, -0.5747595264191646, 0.5334936490538902));
    REQUIRE(+r * +q * q == UnitAxis3D<double>(0.4330127018922193, -0.899519052838329, 0.05801270189221916));
    REQUIRE(+r * r == UnitAxis3D<double>());

    REQUIRE(-m * UnitAxis3D<double>() == m.inv()); 
    REQUIRE(-n * UnitAxis3D<double>() == n.inv()); 
    REQUIRE(-p * UnitAxis3D<double>() == p.inv()); 
    REQUIRE(-q * q == UnitAxis3D<double>());
    REQUIRE(-q * r == UnitAxis3D<double>(0.25, 0.43301270189221935, -0.8660254037844387));
    REQUIRE(-r * q == UnitAxis3D<double>(0.25, -0.058012701892219395, 0.9665063509461096));
    REQUIRE(-r * +r * q == q);
    REQUIRE(-r * -q * q == r.inv());
    REQUIRE(-r * r == UnitAxis3D<double>());
  }

  //SECTION("Boxminus")
  //{
  //  REQUIRE(boxminus(m, m) == Approx(0));
  //  REQUIRE(boxminus(n, n) == Approx(0));
  //  REQUIRE(boxminus(q, r) == Approx(pi/2));
  //  REQUIRE(boxminus(r, q) == Approx(pi/2));
  //  REQUIRE(boxminus(p, q) == Approx(pi/3));
  //  REQUIRE(boxminus(q, p) == Approx(-pi/3));
  //  REQUIRE(boxminus(n, r) == Approx(7*pi/12 - pi));
  //  REQUIRE(boxminus(r, n) == Approx(-7*pi/12 + pi));
  //}

  //SECTION("Boxplus")
  //{
  //  REQUIRE(boxplus(m, 0.0) == m);
  //  REQUIRE(boxplus(m, pi/4) == n);
  //  REQUIRE(boxplus(m, -pi/3) == r);
  //  REQUIRE(boxplus(q, -pi/2) == r);
  //  REQUIRE(boxplus(r, pi/2) == q);
  //}

  //SECTION("Manifold axioms")
  //{
  //  REQUIRE(boxplus(m, 0.0) == m);
  //  REQUIRE(boxplus(r, boxminus(q, r)) == q);
  //  REQUIRE(boxplus(n, boxminus(r, n)) == r);
  //  REQUIRE(boxminus(boxplus(q, 0.55), q) == Approx(0.55));
  //  REQUIRE(boxminus(boxplus(r, -0.32), r) == Approx(-0.32));
  //}

  //SECTION("Distance")
  //{
  //  REQUIRE(distance(m, n) == Approx(pi/4));
  //  REQUIRE(distance(n, m) == Approx(pi/4));
  //  REQUIRE(distance(q, r) == Approx(pi/2));
  //  REQUIRE(distance(r, q) == Approx(pi/2));
  //  REQUIRE(distance(m, m) == Approx(0));
  //  REQUIRE(distance(n, n) == Approx(0));
  //  REQUIRE(distance(q, r) == Approx(pi/2));
  //  REQUIRE(distance(r, q) == Approx(pi/2));
  //  REQUIRE(distance(p, q) == Approx(pi/3));
  //  REQUIRE(distance(q, p) == Approx(pi/3));
  //  REQUIRE(distance(n, r) == Approx(pi - 7*pi/12));
  //  REQUIRE(distance(r, n) == Approx(pi - 7*pi/12));
  //}

  //SECTION("Dot")
  //{
  //  REQUIRE(dot(m, n) == Approx(std::cos(pi/4)));
  //  REQUIRE(dot(n, m) == Approx(std::cos(pi/4)));
  //  REQUIRE(dot(q, r) == Approx(std::cos(pi/2)));
  //  REQUIRE(dot(r, q) == Approx(std::cos(pi/2)));
  //  REQUIRE(dot(m, m) == Approx(std::cos(0)));
  //  REQUIRE(dot(n, n) == Approx(std::cos(0)));
  //  REQUIRE(dot(q, r) == Approx(std::cos(pi/2)));
  //  REQUIRE(dot(r, q) == Approx(std::cos(pi/2)));
  //  REQUIRE(dot(p, q) == Approx(std::cos(pi/3)));
  //  REQUIRE(dot(q, p) == Approx(std::cos(pi/3)));
  //  REQUIRE(dot(n, r) == Approx(std::cos(pi - 7*pi/12)));
  //  REQUIRE(dot(r, n) == Approx(std::cos(pi - 7*pi/12)));
  //}
}

