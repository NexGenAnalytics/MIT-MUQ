#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/ParallelTempering.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;




TEST(MCMC, ParallelTempering_FromOpts) {

    // parameters for the sampler
    pt::ptree pt;
    pt.put("NumSamples", 1e3); // number of Monte Carlo samples
    pt.put("PrintLevel",0);
    pt.put("Inverse Temperatures","0.0,0.25,0.5,0.75,1.0");
    pt.put("Swap Increment", 10);
    pt.put("Swap Type", "DEO");
    pt.put("Adapt Start", 100); // Start adapting after this many steps
    pt.put("Kernel Lists", "Kernel1;Kernel1;Kernel1;Kernel1;Kernel1"); // the transition kernel for the first dimension
    pt.put("Kernel1.Method","MHKernel");
    pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
    pt.put("Kernel1.MyProposal.Method", "MHProposal");
    pt.put("Kernel1.MyProposal.ProposalVariance", 1.0); // the variance of the isotropic MH proposal

    // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
    const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
    auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

    // create a sampling problem
    auto problem = std::make_shared<InferenceProblem>(dist, dist);


    ParallelTempering sampler(pt, problem);

    auto samps = sampler.Run(mu);

    Eigen::VectorXd q = samps->Mean();
    Eigen::VectorXd mcse = samps->StandardError();
    std::cout << "MCSE: " << mcse.transpose() << std::endl;
    EXPECT_NEAR(mu(0), q(0), 3.0*mcse(0));
    EXPECT_NEAR(mu(1), q(1), 3.0*mcse(1));

    q = samps->Variance();
    EXPECT_NEAR(0.5, q(0), 6*mcse(0));
    EXPECT_NEAR(0.5, q(1), 6*mcse(0));
    
}

