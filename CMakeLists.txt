cmake_minimum_required(VERSION 3.3.1)
project(axes)
set(CMAKE_CXX_STANDARD 11) # set standard to c++11

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(include ${EIGEN3_INCLUDE_DIR})

# Finding source files
file(GLOB SOURCES src/*.cpp)
file(GLOB TEST_SOURCES test/*.cpp)

# Add library
add_library(axes INTERFACE)
target_include_directories(axes INTERFACE include)

# Setup Catch for testing
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE test)
add_executable(test_axes ${TEST_SOURCES} ${SOURCES})
target_link_libraries(test_axes Catch)

# Add tests
enable_testing()
add_test(NAME AxesTests COMMAND test_axes)

