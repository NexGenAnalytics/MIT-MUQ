#include <iostream>
#include <string>

#include "HTTPModPiece.h"


int main(int argc, char* argv[])
{
  if (argc <= 1) {
    std::cerr << "Call with: ./binary HOSTNAME:PORT [BEARER_TOKEN]" << std::endl;
    exit(-1);
  }
  std::string host(argv[1]);

  std::string bearer_token = "";
  if (argc > 2) {
    bearer_token = std::string(argv[2]);
  } else {
    std::cout << "No bearer token was passed. Proceeding without token." << std::endl;
  }

  httplib::Headers headers;
  headers.insert(httplib::make_bearer_token_authentication_header(bearer_token));

  HTTPModPiece modPiece(host, headers);


  Eigen::VectorXd input = Eigen::VectorXd::Ones(1);
  input(0) = 0.1;
  std::vector<Eigen::VectorXd> in = {input};

  std::cout << "Evaluating..." << std::endl;
  std::vector<Eigen::VectorXd> out = modPiece.Evaluate(in);
  std::cout << "Received output" << std::endl;
  std::cout << out[0] << std::endl;

  std::cout << "Evaluating..." << std::endl;
  out = modPiece.Evaluate(in);
  std::cout << "Received output" << std::endl;
  std::cout << out[0] << std::endl;

  std::cout << "Request server to shutdown" << std::endl;
  modPiece.ShutdownServer();

  return 0;
}