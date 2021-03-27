#ifndef HTTPMODPIECE
#define HTTPMODPIECE
#include "MUQ/Modeling/ModPiece.h"
#include "HTTPComm.h"

class HTTPModPiece : public muq::Modeling::ModPiece {
public:

  HTTPModPiece(const std::string host)
   : host(host),
     ModPiece(read_input_size(host), read_output_size(host))
   {
     this->outputs.resize(this->numOutputs);
   }

  void ShutdownServer() {
    httplib::Client cli(host.c_str());
    cli.Post("/Quit");
  }

private:

  Eigen::VectorXi read_input_size(const std::string host){
    httplib::Client cli(host.c_str());

    auto res = cli.Get("/GetInputSizes");
    json response_body = json::parse(res->body);
    std::vector<int> outputvec = response_body["inputSizes"].get<std::vector<int>>();
    return stdvector_to_eigenvectori(outputvec);
  }

  Eigen::VectorXi read_output_size(const std::string host){
    httplib::Client cli(host.c_str());

    auto res = cli.Get("/GetOutputSizes");
    json response_body = json::parse(res->body);
    std::vector<int> outputvec = response_body["outputSizes"].get<std::vector<int>>();
    return stdvector_to_eigenvectori(outputvec);
  }

  void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    httplib::Client cli(host.c_str());
    json request_body;

    for (int i = 0; i < this->numInputs; i++) {
      request_body["input" + std::to_string(i)] = eigenvectord_to_stdvector(inputs[i]);
    }
    request_body["level"] = 0;

    auto res = cli.Post("/Evaluate", request_body.dump(), "text/plain");

    for (int i = 0; i < this->numOutputs; i++) {
      json response_body = json::parse(res->body);
      std::vector<double> outputvec = response_body["output" + std::to_string(i)].get<std::vector<double>>();
      outputs[i] = stdvector_to_eigenvectord(outputvec);
    }
  }

  const std::string host;
};

#endif