
#include <iostream>
#include <boost/array.hpp>
#include <boost/asio.hpp>
#include "spdlog/spdlog.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "comm.h"

using boost::asio::ip::tcp;


#include <MUQ/Modeling/ModPiece.h>


class TCPModPiece : public muq::Modeling::ModPiece {
public:

  TCPModPiece(tcp::socket& socket)
   : socket(socket),
     ModPiece(read_input_size(socket), read_output_size(socket))
   {
     this->outputs.resize(this->numOutputs);
   }

private:

  Eigen::VectorXi read_input_size(tcp::socket& socket) {
    send_string(socket, "dimIn");
    return read_vector_i(socket);
  }

  Eigen::VectorXi read_output_size(tcp::socket& socket) {
    send_string(socket, "dimOut");
    return read_vector_i(socket);
  }

  void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    send_string(socket, "sample");
    for (int i = 0; i < this->numInputs; i++)
      send_vector(socket, inputs[i]);
    for (int i = 0; i < this->numOutputs; i++)
      outputs[i] = read_vector(socket);
  }

  tcp::socket& socket;
};


int main(int argc, char* argv[])
{
  try
  {
    if (argc != 2)
    {
      spdlog::error("Usage: client <host>");
      return 1;
    }
    boost::asio::io_service io_service;
    tcp::resolver resolver(io_service);
    tcp::resolver::query query(argv[1], "4242");
    tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
    tcp::resolver::iterator end;
    tcp::socket socket(io_service);
    boost::system::error_code error = boost::asio::error::host_not_found;
    while (error && endpoint_iterator != end)
    {
      socket.close();
      socket.connect(*endpoint_iterator++, error);
    }
    if (error)
      throw boost::system::system_error(error);

    TCPModPiece modPiece(socket);


    const int dim = 4;
    Eigen::VectorXd input = Eigen::VectorXd::Ones(dim);
    input(2) = 42;
    std::vector<Eigen::VectorXd> in = {input};

    std::vector<Eigen::VectorXd> out = modPiece.Evaluate(in);
    out = modPiece.Evaluate(in);
    std::cout << out[0] << std::endl;

    out = modPiece.Evaluate(in);

    for (int i = 0; i < modPiece.numOutputs; i++)
      std::cout << out[i] << std::endl;

    /*
    boost::property_tree::ptree pt = read_ptree(socket);
    std::cout << pt.get<std::string>("command") << std::endl;
    */

    /*for (int i = 0; i < 3; i++) {
      send_string(socket, "sample");
      Eigen::VectorXd myvector = read_vector(socket);
      std::cout << myvector << std::endl;
    }*/
    send_string(socket, "shutdown");

  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
  return 0;
}