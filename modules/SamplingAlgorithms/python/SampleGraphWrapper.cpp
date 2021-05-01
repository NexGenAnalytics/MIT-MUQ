#include "AllClassWrappers.h"

#include <pybind11/pybind11.h>

#include "MUQ/SamplingAlgorithms/SampleGraphs/SampleGraph.h"

#include "MUQ/Utilities/PyDictConversion.h"

namespace py = pybind11;
namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

void muq::SamplingAlgorithms::PythonBindings::SampleGraphWrapper(py::module& m) {
  py::class_<SampleGraph, std::shared_ptr<SampleGraph> > sampleGraph(m, "SampleGraph");
  sampleGraph.def(py::init( [](std::shared_ptr<RandomVariable> const& rv, py::dict const& d) { return new SampleGraph(rv, ConvertDictToPtree(d)); } ));
  sampleGraph.def(py::init( [](Eigen::MatrixXd const& mat, py::dict const& d) { return new SampleGraph(mat, ConvertDictToPtree(d)); } ));
  sampleGraph.def("Point", &SampleGraph::Point);
  sampleGraph.def("SquaredBandwidth", &SampleGraph::SquaredBandwidth);
  sampleGraph.def("NumSamples", &SampleGraph::NumSamples);
  sampleGraph.def("BandwidthParameterCost", &SampleGraph::BandwidthParameterCost);
}
