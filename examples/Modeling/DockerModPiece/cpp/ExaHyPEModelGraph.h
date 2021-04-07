#include <MUQ/Modeling/Distributions/Density.h>
#include <MUQ/Modeling/Distributions/DensityProduct.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>
#include <MUQ/Modeling/OneStepCachePiece.h>
#include <MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h>
#include <MUQ/Modeling/WorkGraph.h>
#include <MUQ/Modeling/IdentityPiece.h>
#include <MUQ/Modeling/ScaleVector.h>

//using namespace muq::Utilities;
using namespace muq::Modeling;
//using namespace muq::SamplingAlgorithms;

std::shared_ptr<WorkGraph> createWorkGraph(std::string host, std::string bearer_token, int level) {

  httplib::Headers headers;
  headers.insert(httplib::make_bearer_token_authentication_header(bearer_token));

  auto model = std::make_shared<HTTPModPiece>(host, headers, Eigen::VectorXi::Ones(1)*2, Eigen::VectorXi::Ones(1)*4, level);

/*
  const int input_dim = 2;
  Eigen::VectorXd prior_mu = Eigen::ArrayXd::Zero(input_dim);
  Eigen::MatrixXd prior_cov = Eigen::MatrixXd::Identity(input_dim, input_dim);
  std::shared_ptr<Density> prior = std::make_shared<Gaussian>(prior_mu, prior_cov)->AsDensity();
*/

  const int output_dim = 4;
  Eigen::VectorXd data = Eigen::ArrayXd::Zero(output_dim);
  data(0) = 1813.9200000008714;
  data(1) = 0.001907;
  data(2) = 5278.9200000007895;
  data(3) = 0.0006368000000000001;
  Eigen::MatrixXd likelihood_cov = Eigen::MatrixXd::Identity(output_dim, output_dim);
  if (level == 0) {
    likelihood_cov(0, 0) = std::pow(300, 2);
    likelihood_cov(1, 1) = std::pow(0.0005, 2);
    likelihood_cov(2, 2) = std::pow(300, 2);
    likelihood_cov(3, 3) = std::pow(0.0005, 2);
  } else if (level == 1) {
    likelihood_cov(0, 0) = std::pow(50, 2);
    likelihood_cov(1, 1) = std::pow(0.0002, 2);
    likelihood_cov(2, 2) = std::pow(50, 2);
    likelihood_cov(3, 3) = std::pow(0.0002, 2);
  } else {
    std::cerr << "WORKGRAPH REQUESTED FOR UNKNOWN LEVEL!" << std::endl;
    exit(-1);
  }
  std::shared_ptr<Density> likelihood = std::make_shared<Gaussian>(data, likelihood_cov)->AsDensity();

  auto parameter = std::make_shared<ScaleVector>(1.0, 2);
  auto qoi = std::make_shared<ScaleVector>(1.0, 2);

  auto graph = std::make_shared<WorkGraph>();

  graph->AddNode(parameter, "parameter");
  graph->AddNode(qoi, "qoi");

  graph->AddNode(model, "model");

  graph->AddNode(likelihood, "likelihood");
  /*graph->AddNode(prior, "prior");
  graph->AddNode(std::make_shared<DensityProduct>(2), "posterior");*/


  //graph->AddEdge("parameter", 0, "prior", 0);
  graph->AddEdge("parameter", 0, "model", 0);
  graph->AddEdge("parameter", 0, "qoi", 0);

  graph->AddEdge("model", 0, "likelihood", 0);

  //graph->AddEdge("prior", 0, "posterior", 0);
  //graph->AddEdge("likelihood", 0, "posterior", 1);

  //graph->Visualize("uq-graph.png");
  return graph;
}