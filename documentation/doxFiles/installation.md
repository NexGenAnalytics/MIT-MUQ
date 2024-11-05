\page muqinstall Installation

<!-- # Getting MUQ
MUQ currently supports Linux and OSX systems.  We do not currently support Windows, but let us know if that is important to you and we may consider it in the future.

On Linux or OSX, there are three primary ways of installing MUQ:
 1. Using the \c muq conda package from the conda-forge channel.  (See instructions below.)
 2. Using the \c muq docker image.
 3. Installing MUQ from source.

# Building from source
 MUQ2 is hosted on <a href=https://bitbucket.org/mituq/muq2>BitBucket</a>.  To get the source code, you can clone our git repository by running
 ```
 git clone https://bitbucket.org/mituq/muq2
 ```
 While it is usually better to clone MUQ so you can easily obtain updates and bugfixes, you can also download a version of MUQ <a href=https://bitbucket.org/mituq/muq2/downloads>here</a>. -->


# Requirements

- MUQ requires a C++17 compiler
- Currently, we only test GNU 11; support or correctness for versions above 11 are currently not guaranteed
- Python 3.12 (The CI runs with the latest python3-dev version for Ubuntu 24.04)
- CMake >= 3.10

# MUQ Groups

## C++ Compile groups

| Compile Group                     | Dependencies                                                        |
|-----------------------------------|---------------------------------------------------------------------|
| `UTILITIES_HDF5`                  | EIGEN3, HDF5, BOOST                                                |
| `UTILITIES_CORE`                  | EIGEN3, BOOST                                                      |
| `UTILITIES_MULTIINDEX`            | EIGEN3, HDF5, BOOST                                                |
| `MODELING_CORE`                   | EIGEN3, BOOST                                                      |
| `MODELING_UMBRIDGE`               | EIGEN3, BOOST                                                      |
| `MODELING_STAN`                   | EIGEN3, BOOST, STANMATH                                            |
| `MODELING_FLANN`                  | EIGEN3, BOOST, NANOFLANN                                           |
| `MODELING_LINEARALGEBRA`          | EIGEN3, BOOST, STANMATH                                            |
| `MODELING_DISTRIBUTIONS`          | EIGEN3, BOOST, STANMATH                                            |
| `MODELING_ODE`                    | EIGEN3, BOOST, STANMATH, SUNDIALS                                  |
| `OPTIMIZATION_CORE`               | EIGEN3, BOOST, STANMATH                                            |
| `OPTIMIZATION_NLOPT`              | EIGEN3, BOOST, STANMATH, NLOPT                                     |
| `APPROXIMATION_REGRESSION`        | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT                    |
| `APPROXIMATION_POLYNOMIALS`       | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                           |
| `APPROXIMATION_QUADRATURE`        | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                           |
| `APPROXIMATION_POLYNOMIALCHAOS`   | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                           |
| `APPROXIMATION_GP_Kernels`        | EIGEN3, BOOST, STANMATH                                            |
| `APPROXIMATION_GP`                | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                           |
| `SAMPLING_PROBLEMS`               | EIGEN3, BOOST, STANMATH                                            |
| `FANCY_SAMPLING_PROBLEMS`         | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT                    |
| `SAMPLING_ALGORITHM`              | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                           |
| `PARALLEL_SAMPLING_ALGORITHM`     | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT, PARCER, SPDLOG, OTF2 |
| `INFERENCE_FILTERING`             | EIGEN3, BOOST, STANMATH                                            |

## Python Compile Groups

| Compile Group                     | Dependencies                                                        |
|-----------------------------------|---------------------------------------------------------------------|
| `UTILITIES_CORE_PYTHON`           | EIGEN3, BOOST, PYTHON                                              |
| `MODELING_CORE_PYTHON`            | EIGEN3, BOOST, PYTHON                                              |
| `MODELING_SUNDIALS_MODELS_PYTHON` | EIGEN3, BOOST, SUNDIALS, PYTHON                                    |
| `OPTIMIZATION_CORE_PYTHON`        | EIGEN3, BOOST, STANMATH, PYTHON                                    |
| `APPROXIMATION_CORE_PYTHON`       | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT, PYTHON            |
| `SAMPLINGALGORITHMS_CORE_PYTHON`  | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT, PYTHON            |
| `INFERENCE_FILTERING_PYTHON`      | EIGEN3, BOOST, STANMATH, PYTHON                                    |



# Dependencies

| Dependency | Version | Download Link                                                                                     |
|------------|---------|---------------------------------------------------------------------------------------------------|
| Boost      | 1.85.0  | [Boost 1.85.0](https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz)               |
| HDF5       | 1.14.1-2| [HDF5 1.14.1-2](https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_1-2/hdf5-1_14_1-2.zip)|
| NLopt      | 2.8.0   | [NLopt 2.8.0](https://github.com/stevengj/nlopt/archive/refs/tags/v2.8.0.zip)                     |
| Sundials   | 5.5.0   | [Sundials 5.5.0](https://github.com/LLNL/sundials/archive/refs/tags/v5.5.0.zip)                   |
| Eigen      | 3.3.7   | [Eigen 3.3.7](https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip)                  |
| nanoflann  | 1.5.5   | [nanoflann 1.5.5](https://github.com/jlblancoc/nanoflann/archive/refs/tags/v1.5.5.zip)            |
| Stan Math  | 2.18.0  | [Stan Math 2.18.0](https://github.com/stan-dev/math/archive/refs/tags/v2.18.0.zip)                |
| Parcer     | commit 3b1ee6dc3d73 | [Parcer](https://bitbucket.org/mituq/parcer/get/3b1ee6dc3d73.zip)                         |
| spdlog     | 1.10.0  | [spdlog 1.10.0](https://github.com/gabime/spdlog/archive/refs/tags/v1.10.0.zip)                   |
| OTF2       | 3.0.3   | [OTF2 3.0.3](https://perftools.pages.jsc.fz-juelich.de/cicd/otf2/tags/otf2-3.0.3/otf2-3.0.3.tar.gz)|

To build the MUQ dependencies, there are two ways:

- You are on your own and use package management, build from source, etc.
- You use the provided build script [link]:

** It is very important to follow the versions defined in the table above; we do not handle version incompatabilities at the moment.**


# Testing MUQ

If MUQ was configured with gtest (e.g., `MUQ_USE_GTEST` was set to `ON`), then compiling MUQ will produce a test executable called `RunAllTests`.   This executable can be run from the `build` directory using
```
./RunAllTests
```
GTest also provides functionality for runnsing a subset of the tests.  The following comman for example, will run all tests with "MCMC" in the name:
```
./RunAllTests --gtest_filter=*MCMC*
```
See the [GoogleTest](http://google.github.io/googletest/advanced.html#running-a-subset-of-the-tests) documentation for more details.

# Configuration Examples

## Default Compile Groups, no MPI

```
cmake -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -S $MUQ_HOME -B $MUQ_BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
    -C $TPL_PFX/tpls_cache.txt \
    -DMUQ_USE_GTEST=ON \
    -DMUQ_ENABLEGROUP_DEFAULT=ON \
    -DMUQ_USE_PYTHON=ON
```

## Default Compile Groups, with MPI

```
cmake -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -S $MUQ_HOME -B $MUQ_BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
    -C $TPL_PFX/tpls_cache.txt \
    -DMUQ_USE_GTEST=ON \
    -DMUQ_USE_MPI=ON \
    -DMUQ_ENABLEGROUP_DEFAULT=ON
```

## Non-default Compile Groups, no MPI

```
cmake -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -S $MUQ_HOME -B $MUQ_BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
    -C $TPL_PFX/tpls_cache.txt \
    -DMUQ_USE_GTEST=ON \
    -DMUQ_ENABLEGROUP_SAMPLING_ALGORITHM=ON \
    -DMUQ_ENABLEGROUP_DEFAULT=OFF \
    -DMUQ_USE_PYTHON=ON
```

# MUQ Python Bindings

To build all Python compile groups:

```
cmake -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -S $MUQ_HOME -B $MUQ_BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
    -C $TPL_PFX/tpls_cache.txt \
    -DMUQ_USE_GTEST=ON \
    -DMUQ_ENABLEGROUP_DEFAULT=ON \
    -DMUQ_USE_PYTHON=ON
```

By default, `DMUQ_USE_PYTHON` is off.

# Using and linking against MUQ

| Compile Group                    | CMake Namespace                  |
|----------------------------------|----------------------------------|
| `UTILITIES_HDF5`                 | `muq::muqUtilities`             |
| `UTILITIES_CORE`                 | `muq::muqUtilities`             |
| `UTILITIES_MULTIINDEX`           | `muq::muqUtilities`             |
| `MODELING_CORE`                  | `muq::muqModeling`              |
| `MODELING_UMBRIDGE`              | `muq::muqModeling`              |
| `MODELING_STAN`                  | `muq::muqModeling`              |
| `MODELING_FLANN`                 | `muq::muqModeling`              |
| `MODELING_LINEARALGEBRA`         | `muq::muqModeling`              |
| `MODELING_DISTRIBUTIONS`         | `muq::muqModeling`              |
| `MODELING_ODE`                   | `muq::muqModeling`              |
| `OPTIMIZATION_CORE`              | `muq::muqOptimization`          |
| `OPTIMIZATION_NLOPT`             | `muq::muqOptimization`          |
| `APPROXIMATION_REGRESSION`       | `muq::muqApproximation`         |
| `APPROXIMATION_POLYNOMIALS`      | `muq::muqApproximation`         |
| `APPROXIMATION_QUADRATURE`       | `muq::muqApproximation`         |
| `APPROXIMATION_POLYNOMIALCHAOS`  | `muq::muqApproximation`         |
| `APPROXIMATION_GP_Kernels`       | `muq::muqApproximation`         |
| `APPROXIMATION_GP`               | `muq::muqApproximation`         |
| `SAMPLING_PROBLEMS`              | `muq::muqSamplingAlgorithms`    |
| `FANCY_SAMPLING_PROBLEMS`        | `muq::muqSamplingAlgorithms`    |
| `SAMPLING_ALGORITHM`             | `muq::muqSamplingAlgorithms`    |
| `PARALLEL_SAMPLING_ALGORITHM`    | `muq::muqSamplingAlgorithms`    |
| `INFERENCE_FILTERING`            | `muq::muqInference`             |


## Example

[source](https://github.com/NexGenAnalytics/MIT-MUQ/blob/master/examples/SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp/CMakeLists.txt)
```
cmake_minimum_required (VERSION 3.10)

project(MetrpolisInGibbs)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(MetropolisInGibbs MetropolisInGibbs.cpp)
target_link_libraries(MetropolisInGibbs muq::muqSamplingAlgorithms)
```