#include "MUQ/Modeling/HTTPModel/HTTPModPiece.h"

using namespace muq::Modeling;


HTTPModPiece::HTTPModPiece(const std::string host, json config, httplib::Headers headers)
: config(config), client(host, headers),
  ModPiece(read_input_size(host, headers), read_output_size(host, headers))
{
  this->outputs.resize(this->numOutputs);
}


Eigen::VectorXi HTTPModPiece::read_input_size(const std::string host, const httplib::Headers& headers){
  // Would prefer to reuse the existing client, circular dependency in constructor though...
  ShallowModPieceClient client(host, headers);
  return client.inputSizes;
}

Eigen::VectorXi HTTPModPiece::read_output_size(const std::string host, const httplib::Headers& headers){
  ShallowModPieceClient client(host, headers);
  return client.outputSizes;
}

void HTTPModPiece::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {
  client.Evaluate(inputs, config);
  outputs = client.outputs;
}
