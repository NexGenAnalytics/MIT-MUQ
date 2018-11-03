#include "MUQ/Approximation/TransportMaps/TransportMap.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Approximation;

TransportMap::TransportMap(unsigned int const totSize) : ModPiece(Eigen::VectorXi::Constant(1, totSize), Eigen::VectorXi::Constant(1, totSize)) {}

void TransportMap::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {
  outputs.resize(1);

  outputs[0] = EvaluateForward(inputs[0]);
}

std::shared_ptr<TransportMap> TransportMap::Identity(unsigned int dim,
                                                     boost::property_tree::ptree options) {

  std::string typeStr = options.get<std::string>("Type");
  TransportMap::IdentityMethodMap const& map = *GetIdentityMethodMap();

  auto iter = map.find(typeStr);
  if(iter != map.end()) {
    return iter->second(dim,options);
  } else {
    std::stringstream msg;
    msg << "TransportMap::Identity: Invalid map type \"" << typeStr << "\".  No \"Identity\" method has been registered for that type.  Available options are:\n";
    for(auto& pair : map)
      msg << "  " << pair.first << "\n";
     throw std::invalid_argument(msg.str());
  }
}

std::shared_ptr<TransportMap> TransportMap::FromSamples(Eigen::MatrixXd const& samps,
                                                        boost::property_tree::ptree options) {

  std::string typeStr = options.get<std::string>("Type");
  TransportMap::SamplesMethodMap const& map = *GetSamplesMethodMap();

  auto iter = map.find(typeStr);
  if(iter != map.end()) {
    return iter->second(samps,options);
  } else {
    std::stringstream msg;
    msg << "TransportMap::FromSamples: Invalid map type \"" << typeStr << "\".  No \"FromSamples\" method has been registered for that type.  Available options are:\n";
    for(auto& pair : map)
      msg << "  " << pair.first << "\n";
    throw std::invalid_argument(msg.str());
  }
}


std::shared_ptr<TransportMap> TransportMap::FromDensity(std::shared_ptr<muq::Modeling::Density> const& dens,
                                                        boost::property_tree::ptree options) {

  std::string typeStr = options.get<std::string>("Type");
  TransportMap::DensityMethodMap const& map = *GetDensityMethodMap();

  auto iter = map.find(typeStr);
  if(iter != map.end()) {
    return iter->second(dens,options);
  } else {
    std::stringstream msg;
    msg << "TransportMap::FromDensity: Invalid map type \"" << typeStr << "\".  No \"Density\" method has been registered for that type.  Available options are:\n";
    for(auto& pair : map)
      msg << "  " << pair.first << "\n";
    throw std::invalid_argument(msg.str());
  }
}

std::shared_ptr<TransportMap::IdentityMethodMap> TransportMap::GetIdentityMethodMap()
{
  static std::shared_ptr<TransportMap::IdentityMethodMap> map;

  if( !map )
    map = std::make_shared<TransportMap::IdentityMethodMap>();

  return map;
}

std::shared_ptr<TransportMap::SamplesMethodMap> TransportMap::GetSamplesMethodMap()
{
  static std::shared_ptr<TransportMap::SamplesMethodMap> map;

  if( !map )
    map = std::make_shared<TransportMap::SamplesMethodMap>();

  return map;
}

std::shared_ptr<TransportMap::DensityMethodMap> TransportMap::GetDensityMethodMap()
{
  static std::shared_ptr<TransportMap::DensityMethodMap> map;

  if( !map )
    map = std::make_shared<TransportMap::DensityMethodMap>();

  return map;
}
