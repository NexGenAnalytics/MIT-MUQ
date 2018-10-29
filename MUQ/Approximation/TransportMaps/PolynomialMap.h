#ifndef POLYNOMIALMAP_H
#define POLYNOMIALMAP_H

#include "MUQ/Approximation/TransportMaps/ConditionableMap.h"

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

#include "MUQ/Modeling/Distributions/Density.h"
#include <boost/property_tree/ptree_fwd.hpp>
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
      PolynomialMap(std::vector<std::shared_ptr<BasisExpansion> > const& expansions,
                    PolynomialMap::InverseMethod const& invMethod = InverseMethod::Comrade);

      virtual ~PolynomialMap() = default;

      virtual Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt, Eigen::VectorXd const& tgtPt0) const override;

      virtual Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const override;

      virtual double LogDeterminant(Eigen::VectorXd const& evalPt) const override;

      //virtual std::shared_ptr<TransportMapBase> Condition(Eigen::VectorXd const& xHead) const override;


      /** Constructs an identity map based on the polynomial basis and orders defined
          in the options ptree.

          <h3>Options:</h3>
          <table>
          <tr><th>Option Key <th> Optional/Required <th> Type <th> Possible Values <th> Default <th> Description
          <tr><td> InverseMethod <td> Optional <td> string <td> Newton, Sturm, Comrade <td> Comrade <td> Specifies the type of solver used to compute polynomial roots when the EvaluateInverse is called.
          <tr><td> IndexSetType <td> Optional <td> string <td> TotalOrder, Hyperbolic <td> TotalOrder <td> Specifies the type of multiindex set to use.
          <tr><td> Order <td> Required <td> int <td> Any nonnegative integer. <td>  -- <td> The maximum polynomial order allowed.  This is used with IndexSetType to define what cross terms are kept in the set.
          <tr><td> NormScale <td> Optional <td> double <td> Any real number in \f$(0,1]\f$. <td> 0.5 <td> Used with the "Hyperbolic" IndexSetType.  For a multiindex \f$\alpha = \left[\alpha_0,\ldots, \alpha_D\right]\f$, the Hyperbolic index set keeps all multiindices satisfying \f$ \left( \sum_{d=1}^D \alpha_d^q\right)^{1/q}\leq a_{max}\f$.  This option specifies \f$q\f$.
          <tr><td> BasisType <td> Required <td> string <td> A name of a child of IndexedScalarBasis <td>  Specifies the type of basis function in the map.  Typical values are "Legendre" and "ProbabilistHermite".
          </table>
      */
      static std::shared_ptr<PolynomialMap> Identity(unsigned int dim,
                                                     boost::property_tree::ptree& options);

      /** Given samples of a target distribution, this function constructs a
          transformation from the target distribution to a standard normal
          reference distribution. Note that the map returned by this function is
          from target->reference, whereas the map constructed by the
          BuildFromDensity function is in the other direction, reference->target.

          @param[in] samps The matrix of samples.  Each column of samps is a sample.

          <table>
          <tr><th>Option Key <th> Optional/Required <th> Type <th> Possible Values <th> Description
          <tr><td> InverseMethod <td> Optional <td> string <td> Newton, Sturm, Comrade <td> Comrade <td> Specifies the type of solver used to compute polynomial roots when the EvaluateInverse is called.
          <tr><td> IndexSetType <td> Optional <td> string <td> TotalOrder, Hyperbolic <td> Specifies the type of multiindex set to use.
          <tr><td> Order <td> Required <td> int <td> Any nonnegative integer. <td> The maximum polynomial order allowed.  This is used with IndexSetType to define what cross terms are kept in the set.
          <tr><td> NormScale <td> Optional <td> double <td> Any real number in \f$(0,1]\f$. <td> Used with the "Hyperbolic" IndexSetType.  For a multiindex \f$\alpha = \left[\alpha_0,\ldots, \alpha_D\right]\f$, the Hyperbolic index set keeps all multiindices satisfying \f$ \left( \sum_{d=1}^D \alpha_d^q\right)^{1/q}\leq a_{max}\f$.  This option specifies \f$q\f$.
          </table>
      */
      static std::shared_ptr<PolynomialMap> FromSamples(Eigen::MatrixXd const& samps,
                                                        boost::property_tree::ptree& options);

      /**
        Given a target density, this function constructs a transformation from a
        standard normal distribution to the target distribution.  Note that the
        map returned by this function is from reference->target, whereas the map
        constructed by the BuildFromSamples function is in the other direction,
        target->reference.
      */
      // static std::shared_ptr<PolynomialMap> FromDensity(std::shared_ptr<muq::Modeling::Density> const& dens,
      //                                                   boost::property_tree::ptree& options);

    protected:

      PolynomialMap(unsigned int dim,
                    boost::property_tree::ptree& options);

      PolynomialMap(Eigen::MatrixXd const& samps,
                    boost::property_tree::ptree& options);

      PolynomialMap(std::shared_ptr<muq::Modeling::Density> const& dens,
                    boost::property_tree::ptree& options);



      void ExtractInverseMethod(boost::property_tree::ptree& options);

      std::vector<std::shared_ptr<muq::Utilities::MultiIndexSet>> ConstructMultis(unsigned int dim,
                                                                  boost::property_tree::ptree& options);

    private:

      void NewtonsMethod(Eigen::VectorXd& result, double const refPt, unsigned int const component) const;

      void PolynomialRootFinding(Eigen::VectorXd& result, double const refPt, unsigned int const component, PolynomialMap::InverseMethod const& invMethod) const;

      std::vector<std::shared_ptr<BasisExpansion> > expansions;

      /// The method we are using to compute the inverse
      InverseMethod invMethod;

      /// Tolerance for Sturm's and Newton's method
      const double tol = 1.0e-12;

      /// Maximum number of iterations for Newton's method
      const unsigned int maxit = 100;
    }; // class PolynomialMap

  }
}

#endif // #ifndef POLYNOMIALMAP_H
