cmake_minimum_required(VERSION 3.10)

# Project name
project(phase-pilot
  VERSION 1.0.0
  DESCRIPTION "Signal detection and tracking software")

# C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_subdirectory(src)

