#ifndef INVERSEMAP_H
#define INVERSEMAP_H

#include "MUQ/Approximation/TransportMaps/TransportMapBase.h"

namespace muq{
  namespace Approximation{

    /** @class InverseMap
        @ingroup TransportMaps
        @brief Defines the inverse of another transport map.
        @details Given a transformation \f$S(x)=r\f$, this class defines the inverse
                 transformation \f$T(r)=x\f$.  It calls the EvaluateInverse function
                 of the original transformation \f$S\f$.
    */
    class InverseMap : public muq::Modeling:TransportMapBase
    {
    public:

      /** Construct the inverse maps from the forward map and an initial condition
          for the inverse solve.
          @param[in] forwardMapIn  The forward transformation \f$S(x)=r\f$.
          @param[in] x0In An initial guess used in the inverse solve.
      */
      InverseMap(std::shared_ptr<TransportMapBase> const& forwardMapIn,
                 Eigen::VectorXd                   const& x0In);

      virtual Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt,
                                              Eigen::VectorXd const& tgtPt0) const override;

      virtual Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const override;

      virtual double LogDeterminant(Eigen::VectorXd const& evalPt) const override;

      virtual std::shared_ptr<TransportMapBase> Inverse(Eigen::VectorXd const& tgtPt0) const override{return forwardMap;};

    private:

      std::shared_ptr<TransportMapBase> forwardMap;
      Eigen::VectorXd x0;

    }; // class InverseTransport

  }
}


#endif // #ifndef INVERSETRANSPORT_H
