#include "MUQ/Approximation/SampleGraphs/SampleGraph.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Approximation;

SampleGraph::SampleGraph(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
samples(SampleRandomVariable(rv, options.get<std::size_t>("NumSamples")))
{
  Initialize(options);
}

SampleGraph::SampleGraph(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, pt::ptree const& options) :
samples(samples)
{
  Initialize(options);
}

void SampleGraph::Initialize(pt::ptree const& options) {
  // get the stride
  const std::size_t stride = NumSamples()/options.get<std::size_t>("Stride", 1);

  // reserve memory
  clouds.reserve(NumSamples()/stride+1);
  kdtrees.reserve(clouds.size());

  // create a bunch of point clouds
  for( std::size_t lag=0; lag<NumSamples(); lag+=stride ) {
    clouds.emplace_back(samples, lag);
    kdtrees.push_back(
      std::make_shared<NanoflannKDTree>(
        StateDimension(),
        *(clouds.end()-1), nanoflann::KDTreeSingleIndexAdaptorParams(options.get<std::size_t>("MaxLeaf", 10))
      ));
  }

  // create the kd trees
  BuildKDTrees();
}

std::shared_ptr<SampleCollection> SampleGraph::SampleRandomVariable(std::shared_ptr<RandomVariable> const& rv, std::size_t const n) {
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }
  return samples;
}

void SampleGraph::BuildKDTrees() const {
  for( const auto& kdtree : kdtrees ) {
    assert(kdtree);
    // (re-)build the kd-tree
    kdtree->buildIndex();
  }
}

Eigen::VectorXd SampleGraph::Point(std::size_t const i) const {
  assert(samples);
  assert(i<samples->size());
  return samples->at(i)->state[0];
}

Eigen::Ref<Eigen::VectorXd> SampleGraph::Point(std::size_t const i) {
  assert(samples);
  assert(i<samples->size());
  return samples->at(i)->state[0];
}

std::size_t SampleGraph::NumSamples() const {
  assert(samples);
  return samples->size();
}

std::size_t SampleGraph::StateDimension() const {
  assert(samples);
  assert(samples->size()>0);
  return samples->at(0)->state[0].size();
}

void SampleGraph::FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, double const radius2, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag) const {
  // make sure the state size matches
  assert(point.size()==StateDimension());

  // unsorted radius set
  nanoflann::RadiusResultSet<double, std::size_t> resultSet(radius2, neighbors); // use squared radius because kd tree metric is the squared euclidean distance

  // choose the kdtree that ignores the first lag samples
  std::size_t ind = clouds.size()-1;
  while( clouds[ind].lag>lag ) { --ind; }
  assert(clouds[ind].kdtree_get_point_count()==kdtrees[ind]->m_size);

  // find the nearest neighbors---neighbors in a specified radius
  kdtrees[ind]->findNeighbors(resultSet, point.data(), nanoflann::SearchParams());

  // add the lag back to the global indcies
  for( auto& neigh : neighbors ) { neigh.first += clouds[ind].lag; }

  // remove the ones that should have been ignored
  auto it = std::remove_if(neighbors.begin(), neighbors.end(), [lag](std::pair<std::size_t, double> const& neigh) { return neigh.first<lag; } );
  neighbors.erase(it, neighbors.end());
}

double SampleGraph::FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag) const {
  // make sure the state size matches
  assert(point.size()==StateDimension());
  assert(k>0);

  // choose the kdtree that ignores the first lag samples
  std::size_t ind = clouds.size()-1;
  while( clouds[ind].lag>lag ) { --ind; }

  // find the nearest neighbors---neighbors in a specified radius
  std::vector<std::size_t> neighInd(k);
  std::vector<double> neighDist(k);
  const std::size_t nfound = kdtrees[ind]->knnSearch(point.data(), k, neighInd.data(), neighDist.data());

  // add the lag back to the global indcies
  neighbors.reserve(k);
  for( std::size_t i=0; i<nfound; ++i ) { neighbors.push_back(std::pair<std::size_t, double>(neighInd[i]+clouds[ind].lag, neighDist[i])); }

  // remove the ones that should have been ignored
  auto it = std::remove_if(neighbors.begin(), neighbors.end(), [lag](std::pair<std::size_t, double> const& neigh) { return neigh.first<lag; } );
  neighbors.erase(it, neighbors.end());

  // compute the average squared distance
  double avg = 0.0;
  int sub = 0 ;
  for( const auto& neigh : neighbors ) {
    if( neigh.second<1.0e-12 ) { ++sub; continue; } // don't include itself
    avg += neigh.second;
  }
  return (sub==neighbors.size()? 0.0 : avg/(neighbors.size()-sub));
}

SampleGraph::Cloud::Cloud(std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> const& samples, std::size_t const lag) :
samples(samples),
lag(lag)
{
  assert(samples->size()>lag);
}

std::size_t SampleGraph::Cloud::kdtree_get_point_count() const {
  assert(samples);
  return samples->size()-lag;
}

double SampleGraph::Cloud::kdtree_get_pt(std::size_t const p, std::size_t const i) const {
  assert(samples);
  return samples->at(lag+p)->state[0] [i];
}
