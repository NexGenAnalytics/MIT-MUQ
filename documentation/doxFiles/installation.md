# Installing MUQ

This guide provides step-by-step instructions for installing MUQ, including setting up required dependencies, configuring compile groups, and testing your installation. Follow these steps carefully to ensure compatibility and successful setup.

---

## System Requirements

To install and use MUQ, you need:

- **Operating System**: Linux or MacOS
- **C++ Compiler**: A C++17-compatible compiler; MUQ is currently tested only with GNU 11. Compatibility with higher versions is not guaranteed.
- **Python**: Python 3.12 (matching the latest `python3-dev` version for Ubuntu 24.04).
- **CMake**: Version 3.10 or higher.

---

## Choosing Compile Groups

MUQ is divided into modular "compile groups" with different dependency requirements. Compile groups allow you to select specific functionalities you want to enable during installation.

If you want to enable the default compile groups (sufficient for most use cases), you can use the `-DMUQ_ENABLEGROUP_DEFAULT=ON` CMake argument when configuring MUQ. This is shown further in the [Configuration Examples section](#configuration-examples).

### C++ Compile Groups

These groups define the primary features and dependencies for C++ compilation.

| Compile Group                     | Dependencies                                                          |
|-----------------------------------|-----------------------------------------------------------------------|
| `UTILITIES_HDF5`                  | EIGEN3, HDF5, BOOST                                                   |
| `UTILITIES_CORE`                  | EIGEN3, BOOST                                                         |
| `UTILITIES_MULTIINDEX`            | EIGEN3, HDF5, BOOST                                                   |
| `MODELING_CORE`                   | EIGEN3, BOOST                                                         |
| `MODELING_UMBRIDGE`               | EIGEN3, BOOST                                                         |
| `MODELING_STAN`                   | EIGEN3, BOOST, STANMATH                                               |
| `MODELING_FLANN`                  | EIGEN3, BOOST, NANOFLANN                                              |
| `MODELING_LINEARALGEBRA`          | EIGEN3, BOOST, STANMATH                                               |
| `MODELING_DISTRIBUTIONS`          | EIGEN3, BOOST, STANMATH                                               |
| `MODELING_ODE`                    | EIGEN3, BOOST, STANMATH, SUNDIALS                                     |
| `OPTIMIZATION_CORE`               | EIGEN3, BOOST, STANMATH                                               |
| `OPTIMIZATION_NLOPT`              | EIGEN3, BOOST, STANMATH, NLOPT                                        |
| `APPROXIMATION_REGRESSION`        | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT                       |
| `APPROXIMATION_POLYNOMIALS`       | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                              |
| `APPROXIMATION_QUADRATURE`        | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                              |
| `APPROXIMATION_POLYNOMIALCHAOS`   | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                              |
| `APPROXIMATION_GP_Kernels`        | EIGEN3, BOOST, STANMATH                                               |
| `APPROXIMATION_GP`                | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                              |
| `SAMPLING_PROBLEMS`               | EIGEN3, BOOST, STANMATH                                               |
| `FANCY_SAMPLING_PROBLEMS`         | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT                       |
| `SAMPLING_ALGORITHM`              | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN                              |
| `PARALLEL_SAMPLING_ALGORITHM`     | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT, PARCER, SPDLOG, OTF2 |
| `INFERENCE_FILTERING`             | EIGEN3, BOOST, STANMATH                                               |

### Python Compile Groups

For building Python bindings, additional dependencies are required.

| Compile Group                     | Dependencies                                            |
|-----------------------------------|---------------------------------------------------------|
| `UTILITIES_CORE_PYTHON`           | EIGEN3, BOOST, PYTHON                                   |
| `MODELING_CORE_PYTHON`            | EIGEN3, BOOST, PYTHON                                   |
| `MODELING_SUNDIALS_MODELS_PYTHON` | EIGEN3, BOOST, SUNDIALS, PYTHON                         |
| `OPTIMIZATION_CORE_PYTHON`        | EIGEN3, BOOST, STANMATH, PYTHON                         |
| `APPROXIMATION_CORE_PYTHON`       | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT, PYTHON |
| `SAMPLINGALGORITHMS_CORE_PYTHON`  | EIGEN3, HDF5, BOOST, STANMATH, NANOFLANN, NLOPT, PYTHON |
| `INFERENCE_FILTERING_PYTHON`      | EIGEN3, BOOST, STANMATH, PYTHON                         |

---

## Dependency Versions

To avoid compatibility issues, follow these exact versions for each dependency:

| Dependency | Version     | Download Link                                                                                      |
|------------|-------------|----------------------------------------------------------------------------------------------------|
| Boost      | 1.85.0      | [Boost 1.85.0](https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz)                |
| HDF5       | 1.14.1-2    | [HDF5 1.14.1-2](https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_1-2/hdf5-1_14_1-2.zip)|
| NLopt      | 2.8.0       | [NLopt 2.8.0](https://github.com/stevengj/nlopt/archive/refs/tags/v2.8.0.zip)                      |
| Sundials   | 5.5.0       | [Sundials 5.5.0](https://github.com/LLNL/sundials/archive/refs/tags/v5.5.0.zip)                    |
| Eigen      | 3.3.7       | [Eigen 3.3.7](https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip)                   |
| nanoflann  | 1.5.5       | [nanoflann 1.5.5](https://github.com/jlblancoc/nanoflann/archive/refs/tags/v1.5.5.zip)             |
| Stan Math  | 2.18.0      | [Stan Math 2.18.0](https://github.com/stan-dev/math/archive/refs/tags/v2.18.0.zip)                 |
| Parcer     | commit 3b1ee6dc3d73 | [Parcer](https://bitbucket.org/mituq/parcer/get/3b1ee6dc3d73.zip)                          |
| spdlog     | 1.10.0      | [spdlog 1.10.0](https://github.com/gabime/spdlog/archive/refs/tags/v1.10.0.zip)                    |
| OTF2       | 3.0.3       | [OTF2 3.0.3](https://perftools.pages.jsc.fz-juelich.de/cicd/otf2/tags/otf2-3.0.3/otf2-3.0.3.tar.gz)|

To install these dependencies, you can either:

- Use package managers or build from source individually, or
- Use the provided build script in [our other repository](https://github.com/NexGenAnalytics/MIT-MUQ-containers): `python build_tpls.py --wdir $PWD --with all`. This will generate a `tpls_cache.txt` file containing all necessary information to give to CMake in order to find the needed dependencies. Using this file with CMake is shown in the next section.

We highly recommend using the build script as MUQ is tested only with this build process.

<blockquote>
**Note**: Ensure you use the exact versions listed to prevent compatibility issues.
</blockquote>

---

## Building and Testing MUQ

### Configuration Examples

Here are some example configurations for building MUQ with and without MPI and Python bindings:

#### Default Compile Groups, No MPI
```bash
cmake -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=Release \
      -S $MUQ_HOME -B $MUQ_BUILD_DIR \
      -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
      -C $TPL_PFX/tpls_cache.txt \
      -DMUQ_USE_GTEST=ON \
      -DMUQ_ENABLEGROUP_DEFAULT=ON \
      -DMUQ_USE_PYTHON=ON
```

#### Default Compile Groups, With MPI
```bash
cmake -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=Release \
      -S $MUQ_HOME -B $MUQ_BUILD_DIR \
      -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
      -C $TPL_PFX/tpls_cache.txt \
      -DMUQ_USE_GTEST=ON \
      -DMUQ_USE_MPI=ON \
      -DMUQ_ENABLEGROUP_DEFAULT=ON
```

#### Non-default Compile Groups, No MPI
```bash
cmake -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=Release \
      -S $MUQ_HOME -B $MUQ_BUILD_DIR \
      -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
      -C $TPL_PFX/tpls_cache.txt \
      -DMUQ_USE_GTEST=ON \
      -DMUQ_ENABLEGROUP_SAMPLING_ALGORITHM=ON \
      -DMUQ_ENABLEGROUP_DEFAULT=OFF \
      -DMUQ_USE_PYTHON=ON
```

### Testing MUQ

If you enabled Google Test (`MUQ_USE_GTEST=ON`), you can run tests with:
```bash
./RunAllTests
```

To run specific tests, use:
```bash
./RunAllTests --gtest_filter=*MCMC*
```
For more options, see the [GoogleTest documentation](http://google.github.io/googletest/advanced.html#running-a-subset-of-the-tests).

---

## MUQ Python Bindings

To build Python bindings, include `-DMUQ_USE_PYTHON=ON`. Here’s an example:
```bash
cmake -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=Release \
      -S $MUQ_HOME -B $MUQ_BUILD_DIR \
      -DCMAKE_INSTALL_PREFIX=$MUQ_INSTALL_DIR \
      -C $TPL_PFX/tpls_cache.txt \
      -DMUQ_USE_GTEST=ON \
      -DMUQ_ENABLEGROUP_DEFAULT=ON \
      -DMUQ_USE_PYTHON=ON
```

<blockquote>
**Note**: Python bindings are disabled by default (`DMUQ_USE_PYTHON=OFF`).
</blockquote>

---

## Linking Against MUQ

Below are some example CMake namespaces for linking with MUQ:

| Compile Group                    | CMake Namespace              |
|----------------------------------|------------------------------|
| `UTILITIES_HDF5`                 | `muq::muqUtilities`          |
| `UTILITIES_CORE`                 | `muq::muqUtilities`          |
| `UTILITIES_MULTIINDEX`           | `muq::muqUtilities`          |
| `MODELING_CORE`                  | `muq::muqModeling`           |
| `MODELING_UMBRIDGE`              | `muq::muqModeling`           |
| `MODELING_STAN`                  | `muq::muqModeling`           |
| `MODELING_FLANN`                 | `muq::muqModeling`           |
| `MODELING_LINEARALGEBRA`         | `muq::muqModeling`           |
| `MODELING_DISTRIBUTIONS`         | `muq::muqModeling`           |
| `MODELING_ODE`                   | `muq::muqModeling`           |
| `OPTIMIZATION_CORE`              | `muq::muqOptimization`       |
| `OPTIMIZATION_NLOPT`             | `muq::muqOptimization`       |
| `APPROXIMATION_REGRESSION`       | `muq::muqApproximation`      |
| `APPROXIMATION_POLYNOMIALS`      | `muq::muqApproximation`      |
| `APPROXIMATION_QUADRATURE`       | `muq::muqApproximation`      |
| `APPROXIMATION_POLYNOMIALCHAOS`  | `muq::muqApproximation`      |
| `APPROXIMATION_GP_Kernels`       | `muq::muqApproximation`      |
| `APPROXIMATION_GP`               | `muq::muqApproximation`      |
| `SAMPLING_PROBLEMS`              | `muq::muqSamplingAlgorithms` |
| `FANCY_SAMPLING_PROBLEMS`        | `muq::muqSamplingAlgorithms` |
| `SAMPLING_ALGORITHM`             | `muq::muqSamplingAlgorithms` |
| `PARALLEL_SAMPLING_ALGORITHM`    | `muq::muqSamplingAlgorithms` |
| `INFERENCE_FILTERING`            | `muq::muqInference`          |

### Example CMake Usage

This example ([source](https://github.com/NexGenAnalytics/MIT-MUQ/blob/master/examples/SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp/CMakeLists.txt)) shows how to link a C++ project with MUQ:

```bash
cmake_minimum_required(VERSION 3.10)
project(MetropolisInGibbs)

set(CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(MetropolisInGibbs MetropolisInGibbs.cpp)
target_link_libraries(MetropolisInGibbs muq::muqSamplingAlgorithms)
```

For more examples, see the [MUQ Examples](https://github.com/NexGenAnalytics/MIT-MUQ/blob/master/examples).

---
