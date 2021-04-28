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
  //sampleGraph.def(py::init( [](std::shared_ptr<SampleCollection> const& samples, py::dict const& d) { return new SampleGraph(samples, ConvertDictToPtree(d)); } ));
  sampleGraph.def("Point", &SampleGraph::Point);
  sampleGraph.def("SquaredBandwidth", &SampleGraph::SquaredBandwidth);
}
