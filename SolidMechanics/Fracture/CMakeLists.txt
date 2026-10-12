cmake_minimum_required(VERSION 3.16)
project(fracture LANGUAGES CXX)

# Only look for these if the parent project hasn't already
if(NOT TARGET mfem)
  find_package(MFEM REQUIRED)
  find_package(MPI REQUIRED)
endif()

# # Path to Eigen
set(EIGEN_inc "/opt/homebrew/opt/eigen/include/eigen3/Eigen")

add_library(fracture
            bilininteg.cpp
            lininteg.cpp
)

add_library(fracture::fracture ALIAS fracture)

target_include_directories(fracture PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN_inc} 
)

target_compile_features(fracture PUBLIC cxx_std_17)

target_link_libraries(fracture PUBLIC
    mfem
    MPI::MPI_CXX
)