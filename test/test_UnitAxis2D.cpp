#include "catch.h"
#include "axis.h"

TEST_CASE("Constructors")
{

  SECTION("Default constructor")
  {
    axes::UnitAxis2D<double> m;
    REQUIRE(m.lambda() == 1.0);
    REQUIRE(m.kappa() == 0.0);
  }

  SECTION("Two parameter constructor")
  {
    axes::UnitAxis2D<double> a(1.0, 0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa() == 0.0);

    axes::UnitAxis2D<double> b(2.0, 0.0);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa() == 0.0);

    axes::UnitAxis2D<double> c(1.0, 1.0);
    REQUIRE(c.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(c.kappa() == Approx(std::sqrt(2.0) / 2.0));

    axes::UnitAxis2D<double> d(-1.0, 1.0);
    REQUIRE(d.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(d.kappa() == Approx(-std::sqrt(2.0) / 2.0));

    axes::UnitAxis2D<double> e(1.0, -1.0);
    REQUIRE(e.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(e.kappa() == Approx(-std::sqrt(2.0) / 2.0));

    axes::UnitAxis2D<double> f(0.0, 0.0);
    REQUIRE(f.lambda() == 1.0);
    REQUIRE(f.kappa() == 0.0);

    axes::UnitAxis2D<double> g(0.0, -1.0);
    REQUIRE(g.lambda() == 0.0);
    REQUIRE(g.kappa() == 1.0);
  }

  SECTION("Angle constructor")
  {
    const long double pi = 3.141592653589793238462643383279502884L;

    axes::UnitAxis2D<double> a(0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa() == 0.0);

    axes::UnitAxis2D<double> b(pi);
    REQUIRE(b.lambda() == Approx(1.0));
    REQUIRE(b.kappa() == Approx(0.0));

    axes::UnitAxis2D<double> c(pi/2.0);
    REQUIRE(c.lambda() == Approx(0.0));
    REQUIRE(c.kappa() == Approx(1.0));

    axes::UnitAxis2D<double> d(pi/3.0);
    REQUIRE(d.lambda() == Approx(0.5));
    REQUIRE(d.kappa() == Approx(std::sqrt(3.0) / 2.0));

    axes::UnitAxis2D<double> e(3.0*pi/4.0);
    REQUIRE(e.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(e.kappa() == Approx(-std::sqrt(2.0) / 2.0));
  }

  SECTION("2-vector constructor")
  {
    Eigen::Vector2d u(1.0, 1.0);
    axes::UnitAxis2D<double> a(u);
    REQUIRE(a.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(a.kappa() == Approx(std::sqrt(2.0) / 2.0));

    Eigen::Vector2d v(0.0, 0.0);
    axes::UnitAxis2D<double> b(v);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa() == 0.0);
  }

  SECTION("Copy constructor")
  {
    axes::UnitAxis2D<double> a(1.0, 1.0);
    axes::UnitAxis2D<double> b(a);
    REQUIRE(a.lambda() == b.lambda());
    REQUIRE(a.kappa() == b.kappa());

    axes::UnitAxis2D<double> c(-1.0, 4.32);
    axes::UnitAxis2D<double> d(c);
    REQUIRE(c.lambda() == d.lambda());
    REQUIRE(c.kappa() == d.kappa());
  }

}

TEST_CASE("Introspection")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  axes::UnitAxis2D<double> m;
  axes::UnitAxis2D<double> n(1.0, 1.0);
  axes::UnitAxis2D<double> p(0.0, -1.0);

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

TEST_CASE("Mathematical methods")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  axes::UnitAxis2D<double> m;
  axes::UnitAxis2D<double> n(1.0, 1.0);
  axes::UnitAxis2D<double> p(0.0, -1.0);
  axes::UnitAxis2D<double> q(2.3, -3.7);

  SECTION("Inv")
  {
    REQUIRE(m.inv() == axes::UnitAxis2D<double>(1.0, 0.0));
    REQUIRE(n.inv() == axes::UnitAxis2D<double>(1.0, -1.0));
    REQUIRE(p.inv() == axes::UnitAxis2D<double>(0.0, 1.0));
    REQUIRE(q.inv() == axes::UnitAxis2D<double>(2.3, 3.7));
  }

  SECTION("Log")
  {
    REQUIRE(m.log() == Approx(0.0));
    REQUIRE(n.log() == Approx(pi/4.0));
    REQUIRE(p.log() == Approx(pi/2.0));
    REQUIRE(q.log() == Approx(-1.014630096674437));
  }
}

TEST_CASE("Operators")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  axes::UnitAxis2D<double> m;
  axes::UnitAxis2D<double> n(1.0, 1.0);
  axes::UnitAxis2D<double> p(0.0, -1.0);

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
    m = axes::UnitAxis2D<double>(4.3, -1.1);
    n = axes::UnitAxis2D<double>(-1.3, -11.9);
    p = axes::UnitAxis2D<double>(7.2, 4.8);

    REQUIRE(m == axes::UnitAxis2D<double>(4.3, -1.1));
    REQUIRE(n == axes::UnitAxis2D<double>(-1.3, -11.9));
    REQUIRE(p == axes::UnitAxis2D<double>(7.2, 4.8));
  }
}

TEST_CASE("Functions")
{
  const long double pi = 3.141592653589793238462643383279502884L;
  axes::UnitAxis2D<double> m;
  axes::UnitAxis2D<double> n(1.0, 1.0);
  axes::UnitAxis2D<double> p(0.0, -1.0);

  SECTION("Equivalency operators")
  {
    axes::UnitAxis2D<double> q;
    axes::UnitAxis2D<double> r(1.0, 1.0);
    axes::UnitAxis2D<double> s(0.0, -1.0);

    REQUIRE(m == m);
    REQUIRE(m == q);
    REQUIRE(n == n);
    REQUIRE(n == r);
    REQUIRE(p == p);
    REQUIRE(p == s);

    REQUIRE(m != n);
    REQUIRE(m != p);
    REQUIRE(m != r);
    REQUIRE(m != s);
  }

  SECTION("Axis product")
  {
  }

  SECTION("Boxminus")
  {
  }

  SECTION("Boxplus")
  {
  }

  SECTION("Distance")
  {
  }

  SECTION("Dot")
  {
  }
}
