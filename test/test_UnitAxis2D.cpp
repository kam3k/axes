#include "catch.h"
#include "axis.h"

TEST_CASE("Constructors")
{

  SECTION("Default constructor")
  {
    axes::UnitAxis2D m;
    REQUIRE(m.lambda() == 1.0);
    REQUIRE(m.kappa() == 0.0);
  }

  SECTION("Two parameter constructor")
  {
    axes::UnitAxis2D a(1.0, 0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa() == 0.0);

    axes::UnitAxis2D b(2.0, 0.0);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa() == 0.0);

    axes::UnitAxis2D c(1.0, 1.0);
    REQUIRE(c.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(c.kappa() == Approx(std::sqrt(2.0) / 2.0));

    axes::UnitAxis2D d(-1.0, 1.0);
    REQUIRE(d.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(d.kappa() == Approx(-std::sqrt(2.0) / 2.0));

    axes::UnitAxis2D e(1.0, -1.0);
    REQUIRE(e.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(e.kappa() == Approx(-std::sqrt(2.0) / 2.0));

    axes::UnitAxis2D f(0.0, 0.0);
    REQUIRE(f.lambda() == 1.0);
    REQUIRE(f.kappa() == 0.0);
  }

  SECTION("Angle constructor")
  {
    const long double pi = 3.141592653589793238462643383279502884L;

    axes::UnitAxis2D a(0.0);
    REQUIRE(a.lambda() == 1.0);
    REQUIRE(a.kappa() == 0.0);

    axes::UnitAxis2D b(pi);
    REQUIRE(b.lambda() == Approx(1.0));
    REQUIRE(b.kappa() == Approx(0.0));

    axes::UnitAxis2D c(pi/2.0);
    REQUIRE(c.lambda() == Approx(0.0));
    REQUIRE(c.kappa() == Approx(1.0));

    axes::UnitAxis2D d(pi/3.0);
    REQUIRE(d.lambda() == Approx(0.5));
    REQUIRE(d.kappa() == Approx(std::sqrt(3.0) / 2.0));

    axes::UnitAxis2D e(3.0*pi/4.0);
    REQUIRE(e.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(e.kappa() == Approx(-std::sqrt(2.0) / 2.0));
  }

  SECTION("2-vector constructor")
  {
    Eigen::Vector2d u(1.0, 1.0);
    axes::UnitAxis2D a(u);
    REQUIRE(a.lambda() == Approx(std::sqrt(2.0) / 2.0));
    REQUIRE(a.kappa() == Approx(std::sqrt(2.0) / 2.0));

    Eigen::Vector2d v(0.0, 0.0);
    axes::UnitAxis2D b(v);
    REQUIRE(b.lambda() == 1.0);
    REQUIRE(b.kappa() == 0.0);
  }

  SECTION("Copy constructor")
  {
    axes::UnitAxis2D a(1.0, 1.0);
    axes::UnitAxis2D b(a);
    REQUIRE(a.lambda() == b.lambda());
    REQUIRE(a.kappa() == b.kappa());

    axes::UnitAxis2D c(-1.0, 4.32);
    axes::UnitAxis2D d(c);
    REQUIRE(c.lambda() == d.lambda());
    REQUIRE(c.kappa() == d.kappa());
  }

}
