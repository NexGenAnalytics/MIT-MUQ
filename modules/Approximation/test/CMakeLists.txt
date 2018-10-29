CreateTestGroup(
  APPROXIMATION_GP_Kernels

  GaussianProcesses/KernelTests.cpp
  #GaussianProcesses/KarhunenLoeveTest.cpp
)

CreateTestGroup(
  APPROXIMATION_GP

  GaussianProcesses/GaussianProcessTests.cpp
  GaussianProcesses/StateSpaceTests.cpp
  GaussianProcesses/KarhunenLoeveTest.cpp
)

CreateTestGroup(
  APPROXIMATION_POLYNOMIALS

  Polynomials/BasisExpansionTest.cpp
  Polynomials/MonotoneExpansionTest.cpp
  Polynomials/PolynomialTests.cpp
  Quadrature/GaussQuadratureTests.cpp
)

CreateTestGroup(
	APPROXIMATION_REGRESSION

	RegressionTests/RegressionTests.cpp
	RegressionTests/LocalRegressionTests.cpp
)

CreateTestGroup(
	APPROXIMATION_TRANSPORT_MAP

	TransportMapTests/PolynomialMapTests.cpp
)

CreateParallelTestGroup(
	APPROXIMATION_REGRESSION

	RegressionTests/LocalRegressionParallelTests.cpp
)
