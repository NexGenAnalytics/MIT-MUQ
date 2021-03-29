#include <iostream>

#include <string>

//#include <resolv.h> // Header included in httplib.h, causing potential issues with Eigen!

#include "HTTPComm.h"

#include <chrono>
#include <thread>


class ExampleModPiece : public ShallowModPiece {
public:

  ExampleModPiece()
   : ShallowModPiece(Eigen::VectorXi::Ones(1)*1, Eigen::VectorXi::Ones(1))
  {
    outputs.push_back(Eigen::VectorXd::Ones(1));
  }

  void Evaluate(std::vector<Eigen::VectorXd> const& inputs) override {
    std::cout << "Entered" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(0));
    const double mu = 0;
    const double sigma = 1;
    outputs[0][0] = - 1.0/2.0 * std::pow(inputs[0][0] - mu, 2) / std::pow(sigma, 2);
    std::cout << "Left" << std::endl;
  }
};

int main(){

  char const* port_cstr = std::getenv("PORT");
  if ( port_cstr == NULL ) {
    std::cerr << "Environment variable PORT not set!" << std::endl;
    exit(-1);
  }


  const int port = atoi(port_cstr);
  ExampleModPiece modPiece;

  serveModPiece(modPiece, "0.0.0.0", port);

  return 0;
}
