add_library(fracture
            bilininteg.cpp
            lininteg.cpp
            energy.cpp
)

add_library(mfemplus::fracture ALIAS fracture)

# Headers and sources live together; the include root is the repository root, so users write
# #include "mfemplus/SolidMechanics/Fracture/bilininteg.hpp". fracture_utils.hpp and solid_mechanics_utils.hpp are
# internal (not included by any public header).
target_include_directories(fracture
    PUBLIC ${PROJECT_SOURCE_DIR}
)

target_compile_features(fracture PUBLIC cxx_std_17)

target_link_libraries(fracture PUBLIC
    mfem
    MPI::MPI_CXX
)
