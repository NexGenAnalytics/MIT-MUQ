#ifndef BANANADISTRIBUTION_H
#define BANANADISTRIBUTION_H

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

namespace muq{
  namespace Modeling{

    /** @class BananaDistribution
        @ingroup Distributions
        @brief Defines a simple banana-shaped distribution in two dimensions.
        @details Let \f$r\sim N(0,1)\f$ be a two dimensional standard normal
        random variable and consider the transformation given by
        \f[
        \left[\begin{array}{c} x_1\\ x_2 \end{array} \right] = \left[\begin{array}{l}a r_1 \\ ba^2(r_1^2 + 1) + r_2/a \end{array}\right].
        \f]
        Combined with the standard normal distribution over \f$r_1,r_2\f$, this transformation induces a distribution over \f$x_1,x_2\f$.
        The BananaDistribution class defines this transformation-induced distribution.
    */
    class BananaDistribution : public Distribution {

    public:

      BananaDistribution(double aIn, double bIn);

      virtual ~BananaDistribution() = default;

    private:

      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual Eigen::VectorXd GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      const double a;
      const double b;

      std::shared_ptr<Distribution> stdNorm;

    }; // class BananaDistribution
  }
}

#endif // # ifndef BANANADISTRIBUTION_H
