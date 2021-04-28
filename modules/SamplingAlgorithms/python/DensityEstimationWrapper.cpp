#include "AllClassWrappers.h"

#include <pybind11/pybind11.h>

#include "MUQ/SamplingAlgorithms/SampleGraphs/DensityEstimation.h"

#include "MUQ/Utilities/PyDictConversion.h"

namespace py = pybind11;
namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

void muq::SamplingAlgorithms::PythonBindings::DensityEstimationWrapper(py::module& m) {
  py::class_<DensityEstimation, SampleGraph, std::shared_ptr<DensityEstimation> > densEstimation(m, "DensityEstimation");
  densEstimation.def(py::init( [](std::shared_ptr<RandomVariable> const& rv, py::dict const& d) { return new DensityEstimation(rv, ConvertDictToPtree(d)); } ));
  //densEstimation.def(py::init( [](std::shared_ptr<SampleCollection> const& samples, py::dict const& d) { return new DensityEstimation(samples, ConvertDictToPtree(d)); } ));
  densEstimation.def("EstimateDensity", &DensityEstimation::EstimateDensity, "Estimate the probability density function.", py::arg("epsilon") = 1.0, py::arg("tune") = true);
}
