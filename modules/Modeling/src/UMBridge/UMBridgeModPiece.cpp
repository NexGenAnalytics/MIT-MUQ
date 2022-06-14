#include "MUQ/Modeling/UMBridge/UMBridgeModPiece.h"

using namespace muq::Modeling;


UMBridgeModPiece::UMBridgeModPiece(const std::string host, json config, httplib::Headers headers)
: config(config), client(host, headers),
  ModPiece(read_input_size(host, headers), read_output_size(host, headers))
{
  this->outputs.resize(this->numOutputs);
}


Eigen::VectorXi UMBridgeModPiece::read_input_size(const std::string host, const httplib::Headers& headers){
  // Would prefer to reuse the existing client, circular dependency in constructor though...
  umbridge::HTTPModel client(host, headers);
  return client.inputSizes;
}

Eigen::VectorXi UMBridgeModPiece::read_output_size(const std::string host, const httplib::Headers& headers){
  umbridge::HTTPModel client(host, headers);
  return client.outputSizes;
}

void UMBridgeModPiece::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {
  if (!client.SupportsEvaluate())
    throw std::runtime_error("Model does not support evaluation!");
  client.Evaluate(inputs, config);
  outputs = client.outputs;
}

void UMBridgeModPiece::GradientImpl(unsigned int outWrt,
                                unsigned int inWrt,
                                muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs,
                                Eigen::VectorXd const& sens) {
  if (client.SupportsGradient()) {
    client.Gradient(outWrt, inWrt, inputs, sens, config);
    gradient = client.gradient;
  } else
    gradient = GradientByFD(outWrt, inWrt, inputs, sens);
}

void UMBridgeModPiece::ApplyJacobianImpl(unsigned int outWrt,
                                    unsigned int inWrt,
                                    muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs,
                                    Eigen::VectorXd const& vec){
  if (client.SupportsApplyJacobian()) {
    client.ApplyJacobian(outWrt, inWrt, inputs, vec, config);
    jacobianAction = client.jacobianAction;
  } else
    jacobianAction = ApplyJacobianByFD(outWrt, inWrt, inputs, vec);
}

void UMBridgeModPiece::ApplyHessianImpl(unsigned int outWrt,
                                    unsigned int inWrt1,
                                    unsigned int inWrt2,
                                    muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs,
                                    Eigen::VectorXd const& sens,
                                    Eigen::VectorXd const& vec){
  if (client.SupportsApplyHessian()) {
    client.ApplyHessian(outWrt, inWrt1, inWrt2, inputs, sens, vec, config);
    hessAction = client.hessAction;
  } else
    hessAction = ApplyHessianByFD(outWrt, inWrt1, inWrt2, inputs, sens, vec);
}