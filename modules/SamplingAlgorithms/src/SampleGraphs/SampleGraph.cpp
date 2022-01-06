#include "MUQ/SamplingAlgorithms/SampleGraphs/SampleGraph.h"

#include <numeric>

#include "MUQ/Optimization/NLoptOptimizer.h"

#include "MUQ/SamplingAlgorithms/SampleGraphs/KernelBandwidthCostFunction.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::SamplingAlgorithms;

SampleGraph::SampleGraph(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
samples(BuildSamples(rv, options.get<std::size_t>("NumSamples"))),
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

SampleGraph::SampleGraph(Eigen::MatrixXd const& mat, pt::ptree const& options) :
samples(BuildSamples(mat)),
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
  ResetIndices();

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

  // set the defaults for the bandwidth parameter optimization
  if( auto opts = options.get_child_optional("BandwidthCostOptimization") ) { bandwidthOptimizationOptions = *opts; }
  bandwidthOptimizationOptions.put("StepSize", bandwidthOptimizationOptions.get<double>("StepSize", 1.0));
  bandwidthOptimizationOptions.put("SparsityTolerance", bandwidthOptimizationOptions.get<double>("SparsityTolerance", 1.0e-1));
  bandwidthOptimizationOptions.put("Ftol.AbsoluteTolerance", bandwidthOptimizationOptions.get<double>("Ftol.AbsoluteTolerance", 1.0e-10));
  bandwidthOptimizationOptions.put("Ftol.RelativeTolerance", bandwidthOptimizationOptions.get<double>("Ftol.RelativeTolerance", 1.0e-10));
  bandwidthOptimizationOptions.put("Xtol.AbsoluteTolerance", bandwidthOptimizationOptions.get<double>("Xtol.AbsoluteTolerance", 1.0e-10));
  bandwidthOptimizationOptions.put("Xtol.RelativeTolerance", bandwidthOptimizationOptions.get<double>("Xtol.RelativeTolerance", 1.0e-10));
  bandwidthOptimizationOptions.put("MaxEvaluations", bandwidthOptimizationOptions.get<std::size_t>("MaxEvaluations", 10000));
  bandwidthOptimizationOptions.put("Algorithm", bandwidthOptimizationOptions.get<std::string>("Algorithm", "SBPLX"));
  bandwidthOptimizationOptions.put("Minimize", false);
}

void SampleGraph::ResetIndices() {
  indices.resize(samples->size());
  std::iota(indices.begin(), indices.end(), 0);
}

std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> SampleGraph::BuildSamples(Eigen::MatrixXd const& mat) {
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<mat.cols(); ++i ) { samples->Add(std::make_shared<SamplingState>(mat.col(i))); }
  return samples;
}

std::shared_ptr<SampleCollection> SampleGraph::BuildSamples(std::shared_ptr<RandomVariable> const& rv, std::size_t const n) {
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
  return samples->at(i)->state[0];
}

Eigen::VectorXd& SampleGraph::Point(std::size_t const i) {
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
  std::vector<unsigned int> neighInd(k);
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
  assert(neighbors.size()==k);

  // return the sum of the squared distance
  double sum = 0.0;
  for( const auto& neigh : neighbors ) { sum += neigh.second; }
  if( sum<1.0e-14 ) {
    //std::cout << "k: " << k << " size: " << neighbors.size() << std::endl;
    //for( const auto& neigh : neighbors ) { std::cout << neigh.second << " "; }
    //std::cout << std::endl << std::endl;
  }
  assert(sum>0.0);
  return sum;
}

Eigen::VectorXd SampleGraph::KernelMatrix(double bandwidthParameter, Eigen::VectorXd const& bandwidth, Eigen::Ref<Eigen::MatrixXd> kernel) const {
  assert(bandwidth.size()==NumSamples());
  assert(kernel.rows()==NumSamples()); assert(kernel.cols()==NumSamples());

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

Eigen::VectorXd SampleGraph::KernelMatrix(double const sparsityTol, double bandwidthParameter, Eigen::VectorXd const& bandwidth, std::vector<Eigen::Triplet<double> >& entries, bool const rebuildTrees) const {
  assert(bandwidthParameter>1.0e-14);
  assert(sparsityTol<1.0+1.0e-14);
  assert(bandwidth.size()==NumSamples());

  // re build the kd trees based on the bandwidth ordering
  if( rebuildTrees ) { BuildKDTrees(bandwidth); }

  // the number of samples
  const std::size_t n = samples->size();

  // compute the entries to the kernel matrix
  entries.clear();
  #pragma omp parallel num_threads(numThreads)
  {
    std::vector<Eigen::Triplet<double> > entries_private;

    //std::cout << "inside parallel" << std::endl;

    assert(indices.size()==n);
    #pragma omp for nowait
    for( std::size_t i=0; i<n; ++i ) {
      //std::cout << std::endl << "i: " << i << " of " << n << std::endl;
      const std::size_t ind = indices[i];

      // we need to find neighbors within this distnace
      double neighDist = bandwidth(ind);
      assert(neighDist>1.0e-14);
      neighDist *= -bandwidthParameter*neighDist*std::log(sparsityTol);
      assert(neighDist>1.0e-14);

      // find the nearest neighbors
      std::vector<std::pair<std::size_t, double> > neighbors;
      FindNeighbors(Point(ind), neighDist, neighbors, i);
      assert(neighbors.size()>0);

      //std::cout << "found neighbors, size: " << neighbors.size() << std::endl;

      // add the kernel evaluations to the kernel matrix
      const double scale = bandwidthParameter*bandwidth(ind);
      //std::size_t cnt = 0;
      std::vector<Eigen::Triplet<double> > entries_neighs;
      entries_neighs.reserve(2*neighbors.size());
      //std::cout << "starting" << std::endl;
      for( const auto& neigh : neighbors ) {
        //std::cout << "cnt: " << cnt++ << " of " << neighbors.size() << std::endl;
	      assert(neigh.first<indices.size());
        const std::size_t jnd = indices[neigh.first];
        assert(i<=neigh.first);
        const double eval = std::exp(-neigh.second/(scale*bandwidth(jnd)));
        //std::cout << "HI" << std::endl;
        if( eval>=sparsityTol ) {
          entries_neighs.emplace_back(ind, jnd, eval);
          if( ind!=jnd ) { entries_neighs.emplace_back(jnd, ind, eval); }
        }
        //std::cout << std::endl;
      }

      //std::cout << "storing entries ... " << std::endl;

      entries_private.insert(entries_private.end(), std::make_move_iterator(entries_neighs.begin()), std::make_move_iterator(entries_neighs.end()));

      //std::cout << "~~~~~~~~~~~~~~~~" << std::endl << std::endl;
    }

    //std::cout << "GOT THIS FAR" << std::endl;

    #pragma omp critical
    entries.insert(entries.end(), std::make_move_iterator(entries_private.begin()), std::make_move_iterator(entries_private.end()));
  }

  //std::cout << "finished computing entries" << std::endl;

  // the sum of each row in the kernel matrix
  Eigen::VectorXd rowsum = Eigen::VectorXd::Zero(n);
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }

  return rowsum;
}


Eigen::VectorXd SampleGraph::KernelMatrix(double const sparsityTol, double bandwidthParameter, Eigen::VectorXd const& bandwidth, Eigen::SparseMatrix<double>& kernel, bool const rebuildTrees) const {
  assert(kernel.rows()==NumSamples()); assert(kernel.cols()==NumSamples());

  // compute the triplets
  std::vector<Eigen::Triplet<double> > entries;
  const Eigen::VectorXd rowsum = KernelMatrix(sparsityTol, bandwidthParameter, bandwidth, entries, rebuildTrees);
  kernel.setFromTriplets(entries.begin(), entries.end());

  return rowsum;
}

double SampleGraph::Kernel(double const scale, Eigen::VectorXd const& xi, Eigen::VectorXd const& xj) const {
  assert(xi.size()==xj.size());
  const Eigen::VectorXd diff = xi-xj;
  return std::exp(-diff.dot(diff)/scale);
}

std::pair<double, double> SampleGraph::TuneKernelBandwidth(Eigen::VectorXd const& bandwidth, double powmin, double powmax, double const epsilon) const {
  BuildKDTrees(bandwidth);

  // create the cost function
  auto cost = std::make_shared<KernelBandwidthCostFunction>(shared_from_this(), bandwidth, bandwidthOptimizationOptions);

  double mincost = cost->BandwidthCost(powmin);
  double maxcost = cost->BandwidthCost(powmax);
  double power = (powmin+powmax)/2.0;
  double bestcost = cost->BandwidthCost(power);

  if( bestcost<mincost ) {
    power = powmin;
    bestcost = mincost;
  } else if( bestcost<maxcost ) {
    power = powmax;
    bestcost = maxcost;
  }

  double prevcost = std::numeric_limits<double>::quiet_NaN();

  while( std::isnan(prevcost) | std::abs(prevcost-bestcost)>1.0e-6 ) {
    prevcost= bestcost;
    //std::cout << "powmax: " << powmax << " powmin: " << powmin << std::endl;
    //std::cout << "mincost: " << mincost << " powmax: " << maxcost << " best cost: " << bestcost << std::endl;

    double lower = (power+mincost)/2.0;
    double lowercost = (std::abs(lower-power)<1.0e-14? bestcost : cost->BandwidthCost(lower));
    if( lowercost>mincost ) {
      powmin = lower;
      mincost = lowercost;
    }
    if( lowercost>bestcost ) {
      power = lower;
      bestcost = lowercost;
    }

    double upper = (power+maxcost)/2.0;
    double uppercost = (std::abs(upper-power)<1.0e-14?  bestcost : cost->BandwidthCost(upper));
    if( uppercost>maxcost ) {
      powmax = upper;
      maxcost = uppercost;
    }
    if( uppercost>bestcost ) {
      power = upper;
      bestcost = uppercost;
    }
  }

  //std::cout << "cost: " << bestcost << " eps: " << std::pow(2.0, power) << std::endl;
  //std::cout << "intrinsic dim: " << 2.0*bestcost << std::endl;
  return std::pair<double, double>(std::pow(2.0, power), bestcost);


  /*const Eigen::VectorXd power = Eigen::VectorXd::LinSpaced(15, powmin, powmax);

  double maxcost = 0.0;
  double optpower = 0.0;
  for( std::size_t i=0; i<power.size(); ++i ) {
    double const cst = cost->BandwidthCost(power(i));
    std::cout << "cost: " << cst << " eps: " << std::pow(2.0, power(i)) << std::endl;
    if( maxcost<cst ) {
      maxcost = cst;
      optpower = power(i);
    }
  }
  std::cout << "intrinsic dim: " << 2.0*maxcost << std::endl;
  return std::pair<double, double>(std::pow(2.0, optpower), maxcost);*/

  auto opt = std::make_shared<NLoptOptimizer>(cost, bandwidthOptimizationOptions);

  // the initial condition for the optimization is the current parameter value
  const std::vector<Eigen::VectorXd> inputs(1, Eigen::VectorXd::Constant(1, std::log2(epsilon)));

  // solve the optimization and update the parameters
  std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);
  const double eps = std::max(1.0e-8, std::pow(2.0, soln.first(0)));
  assert(eps>1.0e-14);
  //std::cout << "cost: " << soln.second << " eps: " << eps << std::endl;
  //std::cout << "intrinsic dim: " << 2.0*soln.second << std::endl;
  return std::pair<double, double>(eps, soln.second);
}

double SampleGraph::BandwidthParameterCost(Eigen::VectorXd const& bandwidth, double const epsilon) const {
  BuildKDTrees(bandwidth);

  auto cost = std::make_shared<KernelBandwidthCostFunction>(shared_from_this(), bandwidth, bandwidthOptimizationOptions);
  return cost->BandwidthCost(log2(epsilon*epsilon));
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

void SampleGraph::WriteToFile(std::string const& filename, std::string const& dataset) const {
  samples->WriteToFile(filename, dataset);
}
