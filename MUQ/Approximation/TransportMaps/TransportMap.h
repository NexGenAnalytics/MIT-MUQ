#ifndef TRANSPORTMAP_H
#define TRANSPORTMAP_H

#include "MUQ/Modeling/ModPiece.h"

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

      TransportMap();

      //virtual Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt,
      //                                        Eigen::VectorXd const& tgtPt0) const = 0;

      //virtual Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const = 0;

      //virtual double LogDeterminant(Eigen::VectorXd const& evalPt) const = 0;


      //virtual std::shared_ptr<TransportMapBase> Inverse(Eigen::VectorXd const& tgtPt0) const;

    private:

      // Implemented in "TransportMap" base class by calling EvaluateForward
      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {}

    }; // class TransportMapBase



  }
}

#endif
