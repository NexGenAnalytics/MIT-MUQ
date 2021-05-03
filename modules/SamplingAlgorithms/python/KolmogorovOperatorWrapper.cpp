#include "AllClassWrappers.h"

#include <pybind11/pybind11.h>

#include "MUQ/SamplingAlgorithms/SampleGraphs/KolmogorovOperator.h"

#include "MUQ/Utilities/PyDictConversion.h"

namespace py = pybind11;
namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

void muq::SamplingAlgorithms::PythonBindings::KolmogorovOperatorWrapper(py::module& m) {
  py::class_<KolmogorovOperator, DensityEstimation, SampleGraph, std::shared_ptr<KolmogorovOperator> > kolmogorov(m, "KolmogorovOperator");
  kolmogorov.def(py::init( [](std::shared_ptr<RandomVariable> const& rv, py::dict const& d) { return new KolmogorovOperator(rv, ConvertDictToPtree(d)); } ));
  kolmogorov.def(py::init( [](Eigen::MatrixXd const& mat, py::dict const& d) { return new KolmogorovOperator(mat, ConvertDictToPtree(d)); } ));
  kolmogorov.def("Eigendecomposition", &KolmogorovOperator::Eigendecomposition, "Estimate the eigendecomposition of the Kolmogorov operator.", py::arg("density"), py::arg("epsilon") = std::numeric_limits<double>::quiet_NaN(), py::arg("tune") = true);
  kolmogorov.def("GradientVectorField", &KolmogorovOperator::GradientVectorField);
  kolmogorov.def("GradientVectorInnerProduct", &KolmogorovOperator::GradientVectorInnerProduct);
  kolmogorov.def_readonly("beta", &KolmogorovOperator::beta);
}
