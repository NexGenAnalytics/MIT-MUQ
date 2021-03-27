#include <iostream>

#include <string>

//#include <resolv.h> // Header included in httplib.h, causing potential issues with Eigen!

#include "HTTPComm.h"



class ExampleModPiece : public ShallowModPiece {
public:

  ExampleModPiece()
   : ShallowModPiece(Eigen::VectorXi::Ones(1)*1, Eigen::VectorXi::Ones(1))
  {
    outputs.push_back(Eigen::VectorXd::Ones(1));
  }

  void Evaluate(std::vector<Eigen::VectorXd> const& inputs) override {
    const double mu = 0;
    const double sigma = 1;
    outputs[0][0] = - 1.0/2.0 * std::pow(inputs[0][0] - mu, 2) / std::pow(sigma, 2);
  }
};

int main(){

  const int port = 4242;
  ExampleModPiece modPiece;

  serveModPiece(modPiece, "0.0.0.0", port);

  return 0;
}
