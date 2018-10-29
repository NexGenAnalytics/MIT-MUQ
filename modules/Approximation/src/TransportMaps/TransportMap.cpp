#include "MUQ/Approximation/TransportMaps/TransportMap.h"

using namespace muq::Approximation;

TransportMap::TransportMap(unsigned int const totSize) : ModPiece(Eigen::VectorXi::Constant(1, totSize), Eigen::VectorXi::Constant(1, totSize)) {}

void TransportMap::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {
  outputs.resize(1);

  outputs[0] = EvaluateForward(inputs[0]);  
}
