#ifndef POLYNOMIALMAP_H
#define POLYNOMIALMAP_H

#include "MUQ/Approximation/TransportMaps/ConditionableMapBase.h"

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"

#include <vector>

namespace muq{
  namespace Approximation{

    /** @class PolynomialMap
        @ingroup TransportMaps
        @brief Defines a lower triangular transport map with multiple polynomial expansions.
    */
    class PolynomialMap : public ConditionableMap
    {
    public:

      virtual Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt,
                                              Eigen::VectorXd const& tgtPt0) const override;

      virtual Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const override;

      virtual double LogDeterminant(Eigen::VectorXd const& evalPt) const override;

      virtual std::shared_ptr<TransportMapBase> Condition(Eigen::VectorXd const& xHead) const override;

    private:

      std::vector<BasisExpansion> expansions;


    }; // class PolynomialMap

  }
}

#endif // #ifndef POLYNOMIALMAP_H
