#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"


#include <iostream>
#include <fstream>

#include <ctime>
#include <string>
#include <boost/asio.hpp>
#include "spdlog/spdlog.h"

using boost::asio::ip::tcp;

using namespace muq::Modeling;


int main(){

  const int dim = 5;
  Eigen::VectorXd input = Eigen::VectorXd::Ones(dim);

  std::string message = "Hello!";

  const int port = 4242;

  try
  {
    boost::asio::io_service io_service;
    tcp::acceptor acceptor(io_service, tcp::endpoint(tcp::v4(), 4242));
    for (;;)
    {
      tcp::socket socket(io_service);
      acceptor.accept(socket);

      boost::system::error_code ignored_error;
      boost::asio::write(socket, boost::asio::buffer(message),
          boost::asio::transfer_all(), ignored_error);

      spdlog::info("Sent");
      break;
    }
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }

  return 0;
}
