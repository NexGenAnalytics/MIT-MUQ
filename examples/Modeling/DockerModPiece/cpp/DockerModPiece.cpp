#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"


#include <iostream>
#include <fstream>

#include <ctime>
#include <string>
#include <boost/asio.hpp>
#include "spdlog/spdlog.h"

#include "comm.h"

using boost::asio::ip::tcp;

using namespace muq::Modeling;


int main(){

  const int dim = 5;
  Eigen::VectorXd dummyoutput = Eigen::VectorXd::Ones(dim);
  dummyoutput(4) = 42;

  const int port = 4242;

  try
  {
    boost::asio::io_service io_service;
    tcp::acceptor acceptor(io_service, tcp::endpoint(tcp::v4(), 4242));
    tcp::socket socket(io_service);
    acceptor.accept(socket);


    while (true) {
      std::string command = read_string(socket);
      if (command == "dimIn") {
        Eigen::VectorXi dimIn = Eigen::VectorXi::Constant(1,4);
        send_vector_i(socket, dimIn);
      } else if (command == "dimOut") {
        Eigen::VectorXi dimOut = Eigen::VectorXi::Constant(2,5);
        send_vector_i(socket, dimOut);
      } else if (command == "sample") {
        spdlog::info("Received sample request");
        Eigen::VectorXd input = read_vector(socket);
        std::cout << "Got input:" << std::endl << input << std::endl;
        send_vector(socket, dummyoutput);
        send_vector(socket, dummyoutput);
        spdlog::info("Sent sample");
      } else if (command == "shutdown") {
        break;
      } else {
        spdlog::error("Received unknown command!");
      }
    }


    /*
    send_string(socket, "{"
    "\"command\":\"Evaluate\",\n"
    "\"content\":\"bla\""
    "}");
    */


    spdlog::info("Quit");
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }

  return 0;
}
