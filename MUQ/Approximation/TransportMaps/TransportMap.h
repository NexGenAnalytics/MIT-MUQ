#ifndef TRANSPORTMAP_H
#define TRANSPORTMAP_H

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include <boost/property_tree/ptree_fwd.hpp>

namespace muq{
  namespace Approximation{

    /** @defgroup TransportMaps
        @ingroup Approximation
        @brief Tools for constructing and using transport maps (i.e., measure transformations)
        @details
    */

    /** @brief Abstract base class for transport maps.
        @ingroup TransportMaps
        @details Consider two random variables \f$x\f$ and \f$r\f$ that both
        take values in \f$\mathbb{R}^N\f$ as well as a nonlinear transport map
        \f$S : \mathbb{R}^N\rightarrow \mathbb{R}^N \f$ that transforms \f$x\f$
        into \f$r\f$, i.e. \f$S(x)=r\f$ in distribution.  This class defines an
        abstract for methods that define \f$S\f$.  There are three functions that
        child classes must implement:
        - EvaluateInverse (Evaluates the inverse map, i.e., \f$x=S^{-1}(r)\f$)
        - EvaluateForward (Evaluates the map, \f$S(x) = r\f$)
        - LogDeterminant (Evaluates the log determinant of the map Jacobian, \f$\det \nabla S (x)\f$)
    */
    class TransportMap : public muq::Modeling::ModPiece {

    public:

      /**
       Constructs an identity map in dim dimensions.  The options enable spaceification
       of the type of map parameterization (e.g., layered, polynomial, monotone)

       <h3>Options:</h3>
       <table>
       <tr><th>Option Key <th> Optional/Required <th> Type <th> Possible Values <th> Description
       <tr><td> Type <td> Required <td> string <td> Polynomial, Monotone, Layered <td> Specifies the type of map parameterization to use.
       <tr><td> OptionsBlock <td> Required <td> string <td> Any string <td> Specifies another block in the ptree defining specific options for the parameterization type.  Options for this block are defined in class corresponding to the specified Type.
       </table>

      */
      static std::shared_ptr<TransportMap> Identity(unsigned int dim,
                                                    boost::property_tree::ptree& options);

      /**
      Constructs a transport map \f$S(x)\f$ that transforms a target random variable
      \f$x\f$ into a standard normal reference random variable \f$r\f$.  In this
      function, only samples of the target distribution are used to construct the map.

      <h3>Options:</h3>
      <table>
      <tr><th>Option Key <th> Optional/Required <th> Type <th> Possible Values <th> Description
      <tr><td> Type <td> Required <td> string <td> Polynomial, Monotone, Layered <td> Specifies the type of map parameterization to use.
      <tr><td> OptionsBlock <td> Required <td> string <td> Any string <td> Specifies another block in the ptree defining specific options for the parameterization type.  Options for this block are defined in class corresponding to the specified Type.
      </table>
      */
      static std::shared_ptr<TransportMap> FromSamples(Eigen::MatrixXd const& samps,
                                                       boost::property_tree::ptree& options);

      /**
        Constructs a transport map \f$T(r)\f$ that transforms a standard normal
        random variable \f$r\f$ into a (non-Gaussian) target density defined by
        the "dens" input variable.  Note that if a gradient-based optimizer is
        specified in the options, then gradients of the target density will be
        needed.

        <h3>Options:</h3>
        <table>
        <tr><th>Option Key <th> Optional/Required <th> Type <th> Possible Values <th> Description
        <tr><td> Type <td> Required <td> string <td> Polynomial, Monotone, Layered <td> Specifies the type of map parameterization to use.
        <tr><td> OptionsBlock <td> Required <td> string <td> Any string <td> Specifies another block in the ptree defining specific options for the parameterization type.  Options for this block are defined in class corresponding to the specified Type.
        </table>
      */
      static std::shared_ptr<TransportMap> FromDensity(std::shared_ptr<muq::Modeling::Density> const& dens,
                                                       boost::property_tree::ptree& options);




      TransportMap(unsigned int const totSize);

      virtual ~TransportMap() = default;

      /**
      @param[in] refPt A point in reference space
      @param[in] tgtPt0 An initial guess for Newton's method
      \return The corresponding point in target space
      */
      virtual Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt, Eigen::VectorXd const& tgtPt0) const = 0;

      virtual Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const = 0;

      virtual double LogDeterminant(Eigen::VectorXd const& evalPt) const = 0;

      //virtual std::shared_ptr<TransportMap> Inverse(Eigen::VectorXd const& tgtPt0) const;

    private:

      // Implemented in "TransportMap" base class by calling EvaluateForward
      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

    }; // class TransportMapBase



  }
}

#endif
