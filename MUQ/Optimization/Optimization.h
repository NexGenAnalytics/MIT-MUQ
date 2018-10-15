#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Optimization/CostFunction.h"

namespace muq {
namespace Optimization {
  /// Solve an optimization problem
  /**
     \f{eqnarray}{
     c &=& \min{J(x; \theta_1, ..., \theta_1)} \        \
     f_i(x) &\leq& 0 \                                  \
     g_i(x) &=& 0
     \f}
  */
  class Optimization : public muq::Modeling::WorkPiece {
  public:

    Optimization(std::shared_ptr<CostFunction> cost,
                 boost::property_tree::ptree const& pt);

    virtual ~Optimization();
    
    /// Add an inequality constraint to the optimization
    /**
       @param[in] ineq The constraint
    */
    virtual void AddInequalityConstraint(std::shared_ptr<CostFunction> ineq)=0;

    /// Add an equality constraint to the optimization
    /**
       NOTE: the NLOPT algorithm used must be able to handle equality constraints
       @param[in] ineq The constraint
    */
    virtual void AddEqualityConstraint(std::shared_ptr<CostFunction> eq)=0;

    /// Solve the optimization problem
    /**
       @param[in] inputs The first input is the variable we are optimizing over, then inputs to the cost function, and inputs to the constraints in the order they were added
       \return First: the argmin, second: the minimum cost
    */
    virtual std::pair<Eigen::VectorXd, double>
    Solve(muq::Modeling::ref_vector<boost::any> const& inputs)=0;

    /// Solve the optimization problem
    /**
       @param[in] args The first input is the variable we are optimizing over, then inputs to the cost function, and inputs to the constraints in the order they were added
       \return First: the argmin, second: the minimum cost
    */
    template<typename ...Args>
    inline std::pair<Eigen::VectorXd, double> Solve(Args... args) {

      Evaluate(args...);

      return std::pair<Eigen::VectorXd, double>(boost::any_cast<Eigen::VectorXd const&>(outputs[0]),
                                                boost::any_cast<double const>(outputs[1]));

    }

  protected:

    /// Evaluate either the cost function or a constraint
    /**
       @param[in] n The size of the input
       @param[in] x The current point
       @param[out] grad The gradient of the cost/constraint 
       @param[in] f_data An Optimization::CostHelper
       \return The cost/constraint value
    */
    static double Cost(unsigned int n,
                       const double* x,
                       double* grad,
                       void* f_data);

    /// Update the inputs if a constraint is added
    /**
       Adding a constraint (potentially) increases the number of inputs to the optimization problem.  If the constraint requires inputs, add them to the optimization.
       @param[in] numNewIns The number of inputs (not the state) that the constraint requires
    */
    virtual void UpdateInputs(unsigned int const numNewIns)=0;
    
    
    /// A structure to help evaluate the cost function and constraints
    struct CostHelper {
      /**
         @param[in] cost The muq::Optimization::CostFunction that evaluates either the cost function or the constraint
         @param[in] firstin The index of the optimziation inputs where this functions inputs begin
      */
      CostHelper(std::shared_ptr<CostFunction> cost, unsigned int const firstin);
      
      virtual ~CostHelper();
      
      /// Given an input to the optimization problem, set the inputs of this fucntion
      /**
         @param[in] ins The inputs to the optimization problem
      */
      void SetInputs(muq::Modeling::ref_vector<boost::any> const& ins);
      
      /// The cost function that we are trying to minimize or a cosntraint
      std::shared_ptr<CostFunction> cost;
      
      /// The index of the optimziation inputs where this functions inputs begin
      const unsigned int firstin;
      
	/// The inputs to this function
      muq::Modeling::ref_vector<Eigen::VectorXd> inputs;
    };

    /// The cost function that we are trying to minimize
    CostHelper opt;

    /// Inequality constraints
    std::vector<CostHelper> ineqConstraints;
    
    /// Equality constraints
    /**
       NOTE: the solver muq::Optimization::Optimization::algorithm must be able to handle equality constraints
    */
    std::vector<CostHelper> eqConstraints;
    
    /// Relative and absoluste tolerances on the cost function value and on the difference between successive values of the state
    const double ftol_rel, ftol_abs, xtol_rel, xtol_abs;
    
    /// Tolerance on the constraints
    const double constraint_tol;
    
    /// Maximum number of cost function evaluations
    const unsigned int maxEvals;
  };

} // namespace Optimization
} // namespace muq

#endif
