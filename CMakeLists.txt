cmake_minimum_required(VERSION 3.15)
# Use SKBUILD_PROJECT_NAME to automatically get the name from pyproject.toml
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)


# Find pybind11, which scikit-build-core ensures is installed
find_package(pybind11 CONFIG REQUIRED)
include(FetchContent)

# 2. Declare the dependency details (name, source, and version tag)
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0 # Specifies the exact version to download
    OVERRIDE_FIND_PACKAGE
)
# Define your compiled extension module. 
# The first argument must be the same name as your project's top-level Python package (splineOPRcpp).
# The path must be relative to the CMakeLists.txt.
pybind11_add_module(splineop_cpp 
    src/py_bindings.cpp
    src/Faulhaber.cpp
    src/QuadraticCost.cpp
    src/SpeedEstimator.cpp
    src/SplineOP.cpp
    src/SplineOP_constrained.cpp
)

# This is the CRITICAL step: it tells CMake to install the compiled file 
# into the Python site-packages directory.
# DESTINATION . ensures it lands in the package root.
install(TARGETS splineop_cpp 
    LIBRARY DESTINATION .
)
FetchContent_MakeAvailable(Eigen)
target_link_libraries(splineop_cpp 
    PRIVATE
        Eigen3::Eigen
)
