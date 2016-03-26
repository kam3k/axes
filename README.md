# Axes
A C++11 implementation of unit axes in both 2D and 3D, tailored for use in state estimation problems. Created by [Marc Gallant](http://kam3k.github.io), originally for use in the [Mining Systems Laboratory](https://msl.engineering.queensu.ca).

## Installation
Axes is a (single) header-only library. As a result, just drop `axes.h` somewhere in your C++ search path. However, it does depend on the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library.

To build the unit tests, Axes uses [Cmake](https://cmake.org). You can view the unit tests in the `tests` directory, which use the popular [Catch](https://github.com/philsquared/Catch) test framework. You do not need to download anything, the Catch header is included in this repository. Make sure the path to your Eigen installation is properly configued in `CMakeLists.txt`, then from the root Axes directory:

```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./test_axes
```

## Background

*Directions* are found in many everyday phenomena. For example, the cardinal directions measured by a compass, the direction of gravitational acceleration, or simply the direction that a car is driving relative to other cars. However, it is sometimes the case that the opposite direction conveys the same information as the direction itself. The quintessential example of this occurrence is the orientation of a plane, such as the (flat) wall of a building. One way to describe the orientation of the wall is its normal vector, which (despite its appearance) is *not* a direction. Because the orientation the wall can equivalently be described by flipping the normal vector to the opposing side of the wall, one must distinguish between the algebraic spaces of normal vectors and directions. The image below shows how an axis **m** is equivalent to –**m** (i.e., flipping **m**) and how it can be used to parameterize the orientation of a line (left) and a plane (right).

![Axes](images/axes.png)

More precisely, an *axis* is an unordered pair of opposing directions. As shown in the above example, a phenomenon that is well-represented by axes is the orientation of planes, which can be parameterized using *three-dimensional* axes. Similarly, the orientation of lines (which are just lower-dimensional planes) can be parameterized using *two-dimensional* axes. This library includes special parameterizations of two-dimensional and three-dimensional axes called *unit axes*, with a focus on using them in state estimation problems. In particular, various operators acting on axes, and the encapsulation of axes as manifolds are explicitly handled.

## Usage

## Issues and Contact

If you discover any problems, I'd appreciate if you submitted them via the issue tracker. Pull requests are certainly welcome! Otherwise, you can find my contact info on my website linked at the beginning of this README.