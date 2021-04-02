#define CPPHTTPLIB_OPENSSL_SUPPORT
#ifndef HTTPMODPIECE
#define HTTPMODPIECE
#include "MUQ/Modeling/ModPiece.h"
#include "HTTPComm.h"

#include <chrono>

class HTTPModPiece : public muq::Modeling::ModPiece {
public:

  HTTPModPiece(const std::string host, httplib::Headers headers, int level = 0)
   : host(host), headers(headers), level(level),
     ModPiece(read_input_size(host, headers), read_output_size(host, headers))
   {
     this->outputs.resize(this->numOutputs);
   }

  HTTPModPiece(const std::string host,
               httplib::Headers headers,
               Eigen::VectorXi const& inputSizes,
               Eigen::VectorXi const& outputSizes,
               int level = 0)
   : host(host), headers(headers), level(level),
     ModPiece(inputSizes, outputSizes)
   {
     this->outputs.resize(this->numOutputs);
   }

  void ShutdownServer() {
    httplib::Client cli(host.c_str());
    cli.Post("/Quit");
  }

private:

  Eigen::VectorXi read_input_size(const std::string host, const httplib::Headers& headers){
    httplib::Client cli(host.c_str());

    std::cout << "GET GetInputSizes" << std::endl;
    auto res = cli.Get("/GetInputSizes", headers);
    std::cout << "got GetInputSizes" << std::endl;
    json response_body = json::parse(res->body);
    std::vector<int> outputvec = response_body["inputSizes"].get<std::vector<int>>();
    return stdvector_to_eigenvectori(outputvec);
  }

  Eigen::VectorXi read_output_size(const std::string host, const httplib::Headers& headers){
    httplib::Client cli(host.c_str());

    std::cout << "GET GetOutputSizes" << std::endl;
    auto res = cli.Get("/GetOutputSizes", headers);
    std::cout << "got GetOutputSizes" << std::endl;
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
    request_body["level"] = level;

    std::cout << "Sending" << std::endl << request_body.dump() << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto res = cli.Post("/Evaluate", headers, request_body.dump(), "text/plain");
    auto current_time = std::chrono::high_resolution_clock::now();

    std::cout << "Got message after " << std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count() << " seconds" << std::endl << res->body << std::endl;

    json response_body = json::parse(res->body);
    for (int i = 0; i < this->numOutputs; i++) {
      std::vector<double> outputvec = response_body["output" + std::to_string(i)].get<std::vector<double>>();
      outputs[i] = stdvector_to_eigenvectord(outputvec);
    }
  }

  int level;
  httplib::Headers headers;
  const std::string host;
};

#endif
