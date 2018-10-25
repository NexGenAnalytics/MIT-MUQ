#ifndef POLYNOMIALMAP_H
#define POLYNOMIALMAP_H

#include "MUQ/Approximation/TransportMaps/ConditionableMap.h"

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

      enum InverseMethod {
        Newton,

        Sturm,

        Comrade
      };

      /**
        @param[in] expansions The basis expensions for each component of the transport map
        @param[in] invMethod The method used to compute the inverse transport map
      */
      PolynomialMap(std::vector<std::shared_ptr<BasisExpansion> > const& expansions, PolynomialMap::InverseMethod const& invMethod = InverseMethod::Comrade);

      virtual ~PolynomialMap() = default;

      virtual Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt, Eigen::VectorXd const& tgtPt0) const override;

      virtual Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const override;

      virtual double LogDeterminant(Eigen::VectorXd const& evalPt) const override;

      //virtual std::shared_ptr<TransportMapBase> Condition(Eigen::VectorXd const& xHead) const override;

    private:

      void NewtonsMethod(Eigen::VectorXd& result, double const refPt, unsigned int const component) const;

      void PolynomialRootFinding(Eigen::VectorXd& result, double const refPt, unsigned int const component, PolynomialMap::InverseMethod const& invMethod) const;

      std::vector<std::shared_ptr<BasisExpansion> > expansions;

      /// The method we are using to compute the inverse
      InverseMethod invMethod;

      /// Tolerance for Newton's method
      const double tol = 1.0e-12;

      /// Maximum number of iterations for Newton's method
      const unsigned int maxit = 100;
    }; // class PolynomialMap

  }
}

#endif // #ifndef POLYNOMIALMAP_H
