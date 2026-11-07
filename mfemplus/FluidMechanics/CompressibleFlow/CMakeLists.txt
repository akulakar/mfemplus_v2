add_library(compressible_flow
            flux.cpp
            limiter.cpp
)

add_library(mfemplus::compressible_flow ALIAS compressible_flow)

# Headers and sources live together; the include root is the repository root, so users write
# #include "mfemplus/FluidMechanics/CompressibleFlow/hyperbolic.hpp"
target_include_directories(compressible_flow
    PUBLIC ${PROJECT_SOURCE_DIR}
)

target_compile_features(compressible_flow PUBLIC cxx_std_17)

target_link_libraries(compressible_flow PUBLIC
    mfem
    MPI::MPI_CXX
)
