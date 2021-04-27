#include "MUQ/Approximation/SampleGraphs/SampleGraph.h"

#include <numeric>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Approximation;

SampleGraph::SampleGraph(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
samples(SampleRandomVariable(rv, options.get<std::size_t>("NumSamples"))),
numThreads(options.get<std::size_t>("NumThreads", 1))
{
  Initialize(options);
}

SampleGraph::SampleGraph(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, pt::ptree const& options) :
samples(samples),
numThreads(options.get<std::size_t>("NumThreads", 1))
{
  Initialize(options);
}

void SampleGraph::Initialize(pt::ptree const& options) {
  // get the stride
  const std::size_t stride = NumSamples()/options.get<std::size_t>("Stride", 1);

  // reserve memory
  clouds.reserve(NumSamples()/stride+1);
  kdtrees.reserve(clouds.size());

  // the initial order is just however they are stored
  indices.resize(samples->size());
  std::iota(indices.begin(), indices.end(), 0);

  // create a bunch of point clouds
  for( std::size_t lag=0; lag<NumSamples(); lag+=stride ) {
    clouds.emplace_back(samples, indices, lag);
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

void SampleGraph::BuildKDTrees(Eigen::VectorXd const& bandwidth) const {
  std::sort(indices.begin(), indices.end(), [&bandwidth](std::size_t const i, std::size_t const j) { return bandwidth(i)>bandwidth(j); });

  BuildKDTrees();
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
  return samples->at(indices[i])->state[0];
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

void SampleGraph::FindNeighbors(Eigen::VectorXd const& point, double const radius2, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag) const {
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

void SampleGraph::FindNeighbors(Eigen::VectorXd const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag) const {
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
}

double SampleGraph::SquaredBandwidth(Eigen::VectorXd const& x, std::size_t const k) const {
  // find the k nearest neighbors
  std::vector<std::pair<std::size_t, double> > neighbors;
  FindNeighbors(x, k, neighbors);

  // return the sum of the squared distance
  double sum = 0.0;
  std::for_each(neighbors.begin(), neighbors.end(), [&sum] (std::pair<std::size_t, double> const& it) { sum += it.second; } );
  return sum;
}

Eigen::VectorXd SampleGraph::KernelMatrix(double bandwidthParameter, Eigen::VectorXd const& bandwidth, Eigen::Ref<Eigen::MatrixXd> kernel) const {
  assert(bandwidth.size()==NumSamples());
  assert(kernel.rows()==NumSamples()); assert(kernel.cols()==NumSamples());

  // compute the squared bandwidth paramter
  bandwidthParameter *= bandwidthParameter;

  #pragma omp parallel for num_threads(numThreads)
  for( std::size_t i=0; i<samples->size(); ++i ) {
    const double scale = bandwidthParameter*bandwidth(i);
    for( std::size_t j=i; j<samples->size(); ++j ) {
      kernel(i, j) = Kernel(scale*bandwidth(j), Point(i), Point(j));
      if( i!=j ) { kernel(j, i) = kernel(i, j); }
    }
  }

  return kernel.rowwise().sum();
}

Eigen::VectorXd SampleGraph::KernelMatrix(double const sparsityTol, double bandwidthParameter, Eigen::VectorXd const& bandwidth, Eigen::SparseMatrix<double>& kernel) const {
  assert(sparsityTol<1.0);
  assert(bandwidth.size()==NumSamples());
  assert(kernel.rows()==NumSamples()); assert(kernel.cols()==NumSamples());

  // compute the squared bandwidth paramter
  bandwidthParameter *= bandwidthParameter;

  // re build the kd trees based on the bandwidth ordering
  BuildKDTrees(bandwidth);

  // the number of samples
  const std::size_t n = samples->size();

  // compute the entries to the kernel matrix
  std::vector<Eigen::Triplet<double> > entries;
  #pragma omp parallel num_threads(numThreads)
  {
    std::vector<Eigen::Triplet<double> > entries_private;

    #pragma omp for nowait
    for( std::size_t i=0; i<n; ++i ) {
      const std::size_t ind = indices[i];

      // we need to find neighbors within this distnace
      double neighDist = bandwidth(ind);
      neighDist *= -bandwidthParameter*neighDist*std::log(sparsityTol);
      assert(neighDist>0.0);

      // find the nearest neighbors
      std::vector<std::pair<std::size_t, double> > neighbors;
      FindNeighbors(Point(ind), neighDist, neighbors, i);

      // add the kernel evaluations to the kernel matrix
      const double scale = bandwidthParameter*bandwidth(ind);
      for( const auto& neigh : neighbors ) {
        const std::size_t jnd = indices[neigh.first];
        assert(i<=neigh.first);
        const double eval = std::exp(-neigh.second/(scale*bandwidth(jnd)));
        if( eval>=sparsityTol ) {
          entries_private.emplace_back(ind, jnd, eval);
          if( ind!=jnd ) { entries_private.emplace_back(jnd, ind, eval); }
        }
      }
    }

    #pragma omp critical
    entries.insert(entries.end(), std::make_move_iterator(entries_private.begin()), std::make_move_iterator(entries_private.end()));
  }
  kernel.setFromTriplets(entries.begin(), entries.end());

  // the sum of each row in the kernel matrix
  Eigen::VectorXd rowsum = Eigen::VectorXd::Zero(n);
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }

  return rowsum;
}

double SampleGraph::Kernel(double const scale, Eigen::VectorXd const& xi, Eigen::VectorXd const& xj) const {
  assert(xi.size()==xj.size());
  const Eigen::VectorXd diff = xi-xj;
  return std::exp(-diff.dot(diff)/scale);
}

SampleGraph::Cloud::Cloud(std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> const& samples, std::vector<std::size_t> const& indices, std::size_t const lag) :
lag(lag),
samples(samples),
indices(indices)
{
  assert(samples->size()>lag);
}

std::size_t SampleGraph::Cloud::kdtree_get_point_count() const {
  assert(samples);
  return samples->size()-lag;
}

double SampleGraph::Cloud::kdtree_get_pt(std::size_t const p, std::size_t const i) const {
  assert(samples);
  return samples->at(indices[lag+p])->state[0] [i];
}
