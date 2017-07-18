#include "MUQ/Modeling/LinearSDE.h"

#include <random>

#include <unsupported/Eigen/MatrixFunctions>

using namespace muq::Modeling;
using namespace muq::Utilities;

LinearSDE::LinearSDE(std::shared_ptr<LinearOperator>    Fin,
                     std::shared_ptr<LinearOperator>    Lin,
                     Eigen::MatrixXd             const& Qin,
                     boost::property_tree::ptree        options) : stateDim(Fin->rows()), F(Fin), L(Lin), Q(Qin), gen(std::random_device()()), normRG(0.0, 1.0)
{
    // Extract options from the ptree
    ExtractOptions(options);
    
    // Compute the Cholesky decomposition of the white noise process
    sqrtQ = Q.llt().matrixL();
};

void LinearSDE::ExtractOptions(boost::property_tree::ptree options)
{
    dt = options.get("SDE.dt",1e-4);
}

Eigen::VectorXd LinearSDE::EvolveState(Eigen::VectorXd const& f0,
                                       double                 T)
{
    Eigen::VectorXd f = f0;

    const int numTimes = std::ceil(T/dt);

    Eigen::VectorXd z(stateDim);
    
    // Take all but the last step.  The last step might be a partial step
    for(int i=0; i<numTimes-1; ++i)
    {
        for(int j=0; j<stateDim; ++j)
            z(j) = normRG(gen);

        z = (sqrt(dt)*sqrtQ*z).eval();
        f += dt*F->Apply(f) + L->Apply( z );
    }

    // Now take the last step
    double lastDt = T-(numTimes-1)*dt;
    for(int j=0; j<stateDim; ++j)
        z(j) = normRG(gen);

    z = (sqrt(lastDt)*sqrtQ*z).eval();
    f += lastDt*F->Apply(f) + L->Apply( z );

    return f;
}


std::pair<Eigen::VectorXd, Eigen::MatrixXd> LinearSDE::EvolveDistribution(Eigen::VectorXd const& mu0,
                                                                          Eigen::MatrixXd const& gamma0,
                                                                          double                 T) const
{

    Eigen::VectorXd mu = mu0;
    Eigen::VectorXd gamma = gamma0;

    const int numTimes = std::ceil(T/dt);

    Eigen::MatrixXd LQLT = L->Apply( L->Apply(Q).transpose().eval() );
    LQLT = 0.5*(LQLT + LQLT.transpose()); // <- Make sure LQLT is symmetric
    
    Eigen::MatrixXd Fgamma;
    
    // Take all but the last step because the last step might be a partial step.
    for(int i=0; i<numTimes-1; ++i)
    {
        mu += dt * F->Apply(mu);
        Fgamma = F->Apply(gamma);
        gamma += dt * (Fgamma + Fgamma.transpose() + LQLT);
    }

    // Take the last step
    double lastDt = T-(numTimes-1)*dt;
    mu += lastDt * F->Apply(mu);
    Fgamma = F->Apply(gamma);
    gamma += lastDt * (Fgamma + Fgamma.transpose() + LQLT);

    return std::make_pair(mu,gamma);
}
