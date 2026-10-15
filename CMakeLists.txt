cmake_minimum_required(VERSION 3.16)
project(mfemplus LANGUAGES CXX)

# Only look for these if the parent project hasn't already
if(NOT TARGET mfem)
  find_package(MFEM REQUIRED)
  find_package(MPI REQUIRED)
endif()

# Headers and sources live together under mfemplus/ and are included as "mfemplus/<module>/...".
add_subdirectory(mfemplus)
