#include <iostream>
#include <string>

#include "HTTPModPiece.h"


int main(int argc, char* argv[])
{

  HTTPModPiece modPiece("localhost:4242");


  Eigen::VectorXd input = Eigen::VectorXd::Ones(1);
  input(0) = 0.1;
  std::vector<Eigen::VectorXd> in = {input};

  std::vector<Eigen::VectorXd> out = modPiece.Evaluate(in);
  out = modPiece.Evaluate(in);
  std::cout << out[0] << std::endl;

  modPiece.ShutdownServer();

  return 0;
}