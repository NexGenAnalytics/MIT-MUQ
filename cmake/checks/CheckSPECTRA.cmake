
# make sure that Eigen supports the "Ref" command
set(CMAKE_REQUIRED_INCLUDES ${SPECTRA_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
    #include <Spectra/SymEigsSolver.h>
    #include <iostream>

    using namespace Spectra;

    int main()
    {
      // We are going to calculate the eigenvalues of M
      Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
      Eigen::MatrixXd M = A + A.transpose();

      // Construct matrix operation object using the wrapper class DenseSymMatProd
      DenseSymMatProd<double> op(M);

      // Construct eigen solver object, requesting the largest three eigenvalues
      SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 3, 6);

      // Initialize and compute
      eigs.init();
      int nconv = eigs.compute();

      // Retrieve results
      Eigen::VectorXd evalues;
      if(eigs.info() == SUCCESSFUL)
          evalues = eigs.eigenvalues();

      return 0;
    }
    "
    SPECTRA_CODE_COMPILES)

set(SPECTRA_TEST_FAIL 1)
set(MUQ_USE_SPECTRA OFF)
if( NOT SPECTRA_CODE_COMPILES )
  set(SPECTRA_TEST_FAIL 0)
  set(MUQ_USE_SPECTRA ON)
endif()
