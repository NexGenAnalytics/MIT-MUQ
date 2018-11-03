#include "AllClassWrappers.h"

#include "MUQ/Approximation/TransportMaps/TransportMap.h"
#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
using namespace muq::Modeling;
using namespace muq::Utilities;

namespace py = pybind11;

void muq::Approximation::PythonBindings::TransportMapWrapper(py::module &m)
{
  py::class_<TransportMap, ModPiece, std::shared_ptr<TransportMap>> tmBase(m, "TransportMap");
  tmBase
    .def_static("Identity", [](unsigned int dim, py::dict d)->std::shared_ptr<TransportMap> {TransportMap::Identity(dim, ConvertDictToPtree(d));})
    .def_static("FromSamples", [](Eigen::MatrixXd const& samps, py::dict d)->std::shared_ptr<TransportMap> {return TransportMap::FromSamples(samps, ConvertDictToPtree(d));})
    .def_static("FromDensity", [](std::shared_ptr<Density> dens, py::dict d)->std::shared_ptr<TransportMap> {return TransportMap::FromDensity(dens, ConvertDictToPtree(d));})
    .def("EvaluateForward", &TransportMap::EvaluateForward)
    .def("EvaluateInverse", &TransportMap::EvaluateInverse)
    .def("LogDeterminant", &TransportMap::LogDeterminant);

}
