#ifndef ODE_H_
#define ODE_H_

#include "MUQ/Modeling/ODEBase.h"

namespace muq {
  namespace Modeling {
    class ODE : public ODEBase {
    public:
      /**
	 The first input is the initial state (at \f$t=0\f$).  It is also the first input to the right hand side muq::Modeling::WorkPiece.  This must either be a double (state size is one) or a vector of doubles.

	 The next set of inputs are the inputs to the right hand side muq::Modeling::WorkPiece.  If the right hand side input takes 2 inputs besides the state, these correspond to inputs 2 and 3 of the muq::Modeling::ODEBase.   Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 Any inputs after the right hand sides inputs are either doubles or vectors of doubles.  The output type is the state at these times.  For example, if the last three inputs (after the RHS inputs) were 2.0 and std::vector<double>({1.0, 2.0}) and 1.5 then there would be three outputs --- the state at time 2.0, a vector of the state at times 1.0 and 2.0, and the state at time 1.5.
	 @param[in] rhs The right hand side of the ODE
	 @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
	 @param[in] algebra A muq::Modeling::AnyAlgebra used to manipulate the state and input parameters (defaults to the MUQ default)
       */
      ODE(std::shared_ptr<WorkPiece> rhs,  boost::property_tree::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra = std::make_shared<AnyAlgebra>());

      virtual ~ODE();

    private:

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      void Integrate(ref_vector<boost::any> const& inputs) const;

      /// Compute the next time to integrate to
      /**
	 @param[out] nextTime first: the next time to integrate to, second: the output index
	 @param[out] timeIndices Each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
	 @param[in] outputTimes We want the state at these times 
       */
      bool NextTime(std::pair<double, int>& nextTime, std::vector<std::pair<unsigned int, unsigned int> >& timeIndices, ref_vector<boost::any> const& outputTimes) const;

      //bool KeepIntegrating(std::vector<unsigned int>& timeIndices, ref_vector<boost::any> const& outputTimes) const;
      
    };
  } // namespace Modeling
} // namespace muq

#endif
