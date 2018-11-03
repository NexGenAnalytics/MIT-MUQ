# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
    APPROXIMATION_CORE_PYTHON
    "Python wrappers for the basic building blocks of the Approximation library."
    "pymuqApproximation"
    "APPROXIMATION_GP;APPROXIMATION_GP_Kernels"
    "PYTHON"
    ""
    Module.cpp
    KernelWrapper.cpp
    GaussianWrapper.cpp
    PolynomialsWrapper.cpp
    KLWrapper.cpp
    QuadratureWrapper.cpp
    TransportMapsWrapper.cpp
)
