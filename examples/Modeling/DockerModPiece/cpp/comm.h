//#include <iostream>

#include <string>
#include <boost/array.hpp>
#include <boost/asio.hpp>
#include "spdlog/spdlog.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Eigen/Core>
using boost::asio::ip::tcp;

void send_string(tcp::socket& socket, const std::string& string) {
    boost::system::error_code ignored_error;
    boost::array<size_t, 1> total_len_buf;
    total_len_buf[0] = string.size();
    socket.write_some(boost::asio::buffer(total_len_buf), ignored_error);
    socket.write_some(boost::asio::buffer(string), ignored_error);

}

void send_vector(tcp::socket& socket, const Eigen::VectorXd& vector) {
    boost::system::error_code ignored_error;
    boost::array<size_t, 1> total_len_buf;
    total_len_buf[0] = vector.size() * sizeof(double);
    socket.write_some(boost::asio::buffer(total_len_buf), ignored_error);
    const double* vc = vector.data();
    socket.write_some(boost::asio::buffer(vc, vector.size() * sizeof(double)), ignored_error);
}

std::string vector_to_string(Eigen::VectorXd& vector) {
  double* vc = vector.data();

  return std::string(reinterpret_cast<char*>(&vc[0]), vector.size() * sizeof(double));
}
Eigen::VectorXd string_to_vector(std::string str) {
  Eigen::VectorXd received = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*>(&str[0]), str.size() * sizeof(char) / sizeof(double));
  return received;
}

void send_vector_i(tcp::socket& socket, const Eigen::VectorXi& vector) {
    boost::system::error_code ignored_error;
    boost::array<size_t, 1> total_len_buf;
    total_len_buf[0] = vector.size() * sizeof(int);
    socket.write_some(boost::asio::buffer(total_len_buf), ignored_error);
    const int* vc = vector.data();
    socket.write_some(boost::asio::buffer(vc, vector.size() * sizeof(int)), ignored_error);
}

// Read full TCP message consisting of a size_t length value and a binary message of that length
std::vector<char> read_buffer(tcp::socket& socket) {
  boost::system::error_code error;
  boost::array<size_t, 1> total_len_buf;
  socket.read_some(boost::asio::buffer(total_len_buf), error);
  size_t total_len = total_len_buf[0];

  std::vector<char> received(total_len);
  std::size_t read_length = 0;
  while (read_length < total_len)
  {
    std::array<char, 128> buf;
    size_t read_size = std::min(size_t(128), total_len - read_length);

    size_t len = socket.read_some(boost::asio::buffer(buf, read_size), error);
    if (error == boost::asio::error::eof) {
      break; // Connection closed cleanly by peer.
    } else if (error)
      throw boost::system::system_error(error); // Some other error.

    for (int i = 0; i < len; i++) {
      received[read_length+i] = buf[i];
    }
    read_length += len;
  }
  return received;
}

std::string read_string(tcp::socket& socket) {
  std::vector<char> buf = read_buffer(socket);

  return std::string(&buf[0], buf.size());
}

Eigen::VectorXd read_vector(tcp::socket& socket) {
  std::vector<char> buf = read_buffer(socket);

  Eigen::VectorXd received = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*>(&buf[0]), buf.size() * sizeof(char) / sizeof(double));
  return received;
}

Eigen::VectorXi read_vector_i(tcp::socket& socket) {
  std::vector<char> buf = read_buffer(socket);

  Eigen::VectorXi received = Eigen::Map<Eigen::VectorXi>(reinterpret_cast<int*>(&buf[0]), buf.size() * sizeof(char) / sizeof(int));
  return received;
}

boost::property_tree::ptree read_ptree(tcp::socket& socket) {
  std::stringstream ss;
  ss << read_string(socket);
  spdlog::info("Got string {}", ss.str());

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  return pt;
}