#include <iostream>

#include <string>

//#include <resolv.h> // Header included in httplib.h, causing potential issues with Eigen!

#include "HTTPComm.h"

#include <chrono>
#include <thread>
#include <iomanip>
#include <stdlib.h>

int test_delay = 0;

class ExampleModPiece : public ShallowModPiece {
public:

  ExampleModPiece()
   : ShallowModPiece(Eigen::VectorXi::Ones(1)*2, Eigen::VectorXi::Ones(1)*4)
  {
    outputs.push_back(Eigen::VectorXd::Ones(4));
  }

  void Evaluate(std::vector<Eigen::VectorXd> const& inputs, json config) override {
    int level = config.value<int>("level", 0);

    std::cout << "Entered for level " << level << std::endl;

    std::ofstream inputsfile ("/tmp/inputs.txt");
    typedef std::numeric_limits<double> dl;
    inputsfile << std::fixed << std::setprecision(dl::digits10);
    for (int i = 0; i < inputs[0].rows(); i++) {
      inputsfile << inputs[0](i) << std::endl;
    }
    inputsfile.close();

    // TODO: Need level dependence!
    int status;
    if(level == 0) {
      status = system("cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l0 && ./ExaHyPE-SWE ../SWE_asagi_limited_l0.exahype2");
    } else if(level == 1) {
      status = system("cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l1 && ./ExaHyPE-SWE ../SWE_asagi_limited_l1.exahype2");
    } else if(level == 2) {
      status = system("cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l2 && ./ExaHyPE-SWE ../SWE_asagi_limited_l2.exahype2");
    } else {
      std::cerr << "Unknown model requested by client!" << std::endl;
      exit(-1);
    }
    std::cout << "Exahype exit status " << status << std::endl;


    std::ifstream outputsfile("/tmp/outputs.txt");
    for (int i = 0; i < outputs[0].rows(); i++) {
      outputsfile >> outputs[0](i);
    }
    outputsfile.close();
    std::cout << "Read outputs from exahype:" << outputs[0] << std::endl;

    std::cout << "Left" << std::endl;
  }
};

int main(){

  char const* port_cstr = std::getenv("PORT");
  if ( port_cstr == NULL ) {
    std::cerr << "Environment variable PORT not set!" << std::endl;
    exit(-1);
  }

  char const* delay_cstr = std::getenv("TEST_DELAY");
  if ( delay_cstr != NULL ) {
    test_delay = atoi(delay_cstr);
  }


  const int port = atoi(port_cstr);
  ExampleModPiece modPiece;

  serveModPiece(modPiece, "0.0.0.0", port);

  return 0;
}
