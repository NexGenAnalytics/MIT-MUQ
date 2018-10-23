#ifndef CONDITIONABLEMAP_H
#define CONDITIONABLEMAP_H

#include "MUQ/Approximation/TransportMaps/TransportMapBase.h"

namespace muq{
  namespace Approximation{

    /** @class ConditionableMap
        @ingroup TransportMaps
        @brief Defines transport maps that can be conditioned on part of the input.
        @details Consider a target random variable \f$x\f$ that is defined through
                 two blocks \f$x= [x_1,x_2]\f$ and a corresponding transport map
                 that is also defined blockwise:
                 \f[
                 \left[\begin{array}{l} S_1(x_1)\\ S_2(x_1,_2)\end{array}\right]  = \left[ \begin{array}{c} r_1\\ r_2\end{array}\right].
                 \f]
                 If \f$x_1\f$ was known, it is then possible to "Condition" the
                map on \f$x_1\f$ and obtain a new map that transforms \f$x_c = x_2 | x_1\f$
                into \f$r_2\f$, i.e.,
                \f[
                \bar{S}(x_c) = S_2(x_1, x_c) = r_2.
                \f]
                The ConditionableMap class defines transformations, typically
                lower triangular, where this type of conditional transformation
                is possible.  In particular, this class adds a "Condition"
                function that accepts the value of \f$x_1\f$ and returns the
                map \f$\bar{S}\f$.
    */
    class ConditionableMap: public muq::Modeling:TransportMapBase
    {
    public:

      virtual std::shared_ptr<TransportMapBase> Condition(Eigen::VectorXd const& xHead) const = 0;


    }; // class ConditionableMap

  }
}

#endif // #ifndef CONDITIONALMAP_H
