#ifndef MAPFACTORY_H
#define MAPFACTORY_H

#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

#include "MUQ/Modeling/Distributions/Density.h"

#include <boost/property_tree/ptree_fwd.hpp>


namespace muq {
  namespace Approximation {

    /** @class MapFactory
        @ingroup TransportMaps
        @brief Factory class with methods for building transport maps.

        <h3>Polynomial Map Options:</h3>
        <table>
        <tr><th>Option Key <th> Optional/Required <th> Type <th> Default Value <th> Possible Values <th> Description
        <tr><td> Order <td> Required <td> int <td> NA <td> Any natural number <td> In conjunction with the LimitType option, this specifies the order of the polynomials to include in the expansions.
        <tr><td> LimitType <td> Optional <td> string <td> "TotalOrder" <td> TotalOrder, MaxOrder, Hyperbolic <td> Specifies the type of limit imposed on the polynomials in the expansion.
        </table>

        <h3>Monotone Map Options:</h3>
        \f[
        S_{i}(x) = g_{i}(x_{1},\ldots, x_{i-1}) + \int_0^{x_i} \left( h_i(x_1,\ldots,x_{i-1}, z) \right) dz
        \f]
        <table>
        <tr><th>Option Key <th> Optional/Required <th> Type <th> Default Value <th> Possible Values <th> Description
        <tr><td> PolyOrder <td> Required <td> int <td> NA <td> Any natural number <td> Defines the polynomial terms used in the expansion for \f$g_i(x_1,\ldots, x_{i-1})\f$.
        <tr><td> PolyLimitType <td> Optional <td> string <td> "TotalOrder" <td> TotalOrder, MaxOrder, Hyperbolic <td> Specifies the type of limit imposed on the polynomials in \f$g_i(x_1,\ldots, x_{i-1})\f$.
        <tr><td> MonoOrder <td> Required <td> int <td> NA <td> Any natural number <td> Defines the polynomial terms used in the expansion for \f$h_i(x_1,\ldots, x_{i-1},z)\f$
        <tr><td> MonoLimitType <td> Optional <td> string <td> "TotalOrder" <td> TotalOrder, MaxOrder, Hyperbolic <td> Specifies the type of limit imposed on the polynomial terms in \f$h_i(x_1,\ldots, x_{i-1},z)\f$.
        </table>

        <h3>Layered Map Options:</h3>
        Something will go here when we implement layered maps.
    */
    class MapFactory {

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
      */
      //static std::shared_ptr<TransportMap> FromDensity(std::shared_ptr<muq::Modeling::Density> const& dens,
      //                                                 boost::property_tree::ptree& options);




    }; // class MapFactory

  }
}

#endif // #ifndef MAPFACTORY_H
