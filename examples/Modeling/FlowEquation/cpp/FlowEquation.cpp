#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/LinearAlgebra/HessianOperator.h"
#include "MUQ/Modeling/LinearAlgebra/StochasticEigenSolver.h"

#include <boost/property_tree/ptree.hpp>

class FlowEquation : public muq::Modeling::ModPiece
{
public:

  FlowEquation(Eigen::VectorXd const& sourceTerm) : muq::Modeling::ModPiece({int(sourceTerm.size())},
                                                                            {int(sourceTerm.size()+1)}) 
  {
    numCells = sourceTerm.size();
    numNodes = sourceTerm.size()+1;

    xs = Eigen::VectorXd::LinSpaced(numNodes,0,1);
    dx = xs(1)-xs(0);

    rhs = BuildRhs(sourceTerm);
  };

protected:

  void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override
  {
    // Extract the conductivity vector from the inptus
    auto& cond = inputs.at(0).get();

    // Build the stiffness matrix
    auto K = BuildStiffness(cond);

    // Solve the sparse linear system and store the solution in the outputs vector
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    
    outputs.resize(1);
    outputs.at(0) = solver.solve(rhs);
  };

  virtual void GradientImpl(unsigned int outWrt,
                            unsigned int inWrt,
                            muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs,
                            Eigen::VectorXd const& sens) override
  {
    // Extract the conductivity vector from the inptus
    auto& cond = inputs.at(0).get();

    // Build the stiffness matrix
    auto K = BuildStiffness(cond);

    // Factor the stiffness matrix for forward and adjoint solves 
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);

    // Solve the forward problem
    Eigen::VectorXd sol = solver.solve(rhs);

    // Solve the adjoint problem
    Eigen::VectorXd adjRhs = BuildAdjointRHS(sens);
    Eigen::VectorXd adjSol = solver.solve(adjRhs);

    // Compute the gradient from the adjoint solution
    gradient = (1.0/dx)*(sol.tail(numCells)-sol.head(numCells)).array() * (adjSol.tail(numCells) - adjSol.head(numCells)).array();
  }

  virtual void ApplyHessianImpl(unsigned int outWrt,
                                 unsigned int inWrt1,
                                 unsigned int inWrt2,
                                 muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs,
                                 Eigen::VectorXd const& sens,
                                 Eigen::VectorXd const& vec) override
  {
    // Extract the conductivity vector from the inptus
    auto& cond = inputs.at(0).get();

    // Build the stiffness matrix
    auto K = BuildStiffness(cond);

    // Factor the stiffness matrix for forward and adjoint solves 
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);

    // Solve the forward problem
    Eigen::VectorXd sol = solver.solve(rhs);

    // If we're using the Hessian $\nabla_{kk} J$
    if((inWrt1==0)&&(inWrt2==0)){

      // Solve the adjoint problem
      Eigen::VectorXd adjRhs = BuildAdjointRHS(sens);
      Eigen::VectorXd adjSol = solver.solve(adjRhs);

      // Solve the incremental forward problem
      Eigen::VectorXd incrForRhs = BuildIncrementalRhs(sol, vec);
      Eigen::VectorXd incrForSol = solver.solve(incrForRhs);

      // Solve the incremental adjoint problem 
      Eigen::VectorXd incrAdjRhs = BuildIncrementalRhs(adjSol, vec);
      Eigen::VectorXd incrAdjSol = solver.solve(incrAdjRhs);

      // Construct the Hessian action
      auto solDeriv = (sol.tail(numCells)-sol.head(numCells))/dx;
      auto adjDeriv = (adjSol.tail(numCells)-adjSol.head(numCells))/dx;
      auto incrForDeriv = (incrForSol.tail(numCells) - incrForSol.head(numCells))/dx;
      auto incrAdjDeriv = (incrAdjSol.tail(numCells) - incrAdjSol.head(numCells))/dx;
      
      hessAction = -(incrAdjDeriv.array() * solDeriv.array() + incrForDeriv.array() * adjDeriv.array());

    // If we're using the mixed Hessian $\nabla_{ks} J$
    }else if((inWrt1==0)&&(inWrt2==1)){

      Eigen::VectorXd temp = solver.solve(vec);
      auto solDeriv = (sol.tail(numCells) - sol.head(numCells))/dx;
      auto tempDeriv = (temp.tail(numCells)-temp.head(numCells))/dx;
      
      hessAction = -dx * solDeriv.array() * tempDeriv.array();

    // We should never see any other options...
    }else{
      assert(false);
    }
  }

  /** Construct the right hand side of the forward problem given a vector containing the source term f_i in each grid cell. */
  Eigen::VectorXd BuildRhs(Eigen::VectorXd const& sourceTerm) const
  {
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(numNodes);

    rhs.segment(1,numNodes-2) = 0.5*dx*(sourceTerm.tail(numNodes-2) + sourceTerm.head(numNodes-2));
    rhs(numNodes-1) = 0.5*dx*sourceTerm(numNodes-2);
    
    return rhs;
  }

  /** Construct the right hand side vector for the adjoint problem. */
  Eigen::VectorXd BuildAdjointRHS(Eigen::VectorXd const& sensitivity) const
  {
    Eigen::VectorXd rhs = -1.0*sensitivity;
    rhs(0) = 0.0; // <- To enforce Dirichlet BC 
    return rhs;
  }

  /** Constructs the right hand side vector for both the incremental forward and incremental adjoint problems. */
  Eigen::VectorXd BuildIncrementalRhs(Eigen::VectorXd const& sol, Eigen::VectorXd const& khat)
  {
    // Compute the derivative of the solution in each cell
    Eigen::VectorXd solGrad = (sol.tail(numCells)-sol.head(numCells))/dx;
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(numNodes);

    unsigned int leftNode, rightNode;
    for(unsigned int cellInd=0; cellInd<numCells; ++cellInd)
    {
      leftNode = cellInd;
      rightNode = cellInd + 1;

      rhs(leftNode) -= dx*khat(cellInd)*solGrad(cellInd);
      rhs(rightNode) += dx*khat(cellInd)*solGrad(cellInd);
    }

    rhs(0) = 0.0; // <- To enforce Dirichlet BC at x=0
    return rhs;
  }

  /** Build the sparse stiffness matrix given a vector of conductivities in each cell. */
  Eigen::SparseMatrix<double> BuildStiffness(Eigen::VectorXd const& condVals) const{

    typedef Eigen::Triplet<double> T;
    std::vector<T> nzVals;

    // Add a large number to K[0,0] to enforce Dirichlet BC
    nzVals.push_back( T(0,0,1e10) );

    unsigned int leftNode, rightNode;
    for(unsigned int cellInd=0; cellInd<numCells; ++cellInd)
    {
      leftNode  = cellInd;
      rightNode = cellInd+1;

      nzVals.push_back( T(leftNode,  rightNode, -condVals(cellInd)/dx) );
      nzVals.push_back( T(rightNode, leftNode,  -condVals(cellInd)/dx) );
      nzVals.push_back( T(rightNode, rightNode,  condVals(cellInd)/dx) );
      nzVals.push_back( T(leftNode,  leftNode,   condVals(cellInd)/dx) );
    }

    Eigen::SparseMatrix<double> stiffMat(numNodes,numNodes);
    stiffMat.setFromTriplets(nzVals.begin(), nzVals.end());

    return stiffMat;
  }


private:

  // Store "mesh" information.  xs contains the node locations and dx is the uniform spacing between nodes
  Eigen::VectorXd xs;
  double dx;
  unsigned int numCells;
  unsigned int numNodes;

  // Will store the precomputed RHS for the forward problem
  Eigen::VectorXd rhs;

}; // end of class SimpleModel



int main(){

  using namespace muq::Modeling;
  using namespace muq::Utilities;

  unsigned int numCells = 20;

  // Vector containing the location of each node in the "mesh"
  Eigen::VectorXd nodeLocs = Eigen::VectorXd::LinSpaced(numCells+1,0,1);

  // Vector containing the midpoint of each cell
  Eigen::VectorXd cellLocs = 0.5*(nodeLocs.head(numCells) + nodeLocs.tail(numCells));

  // The source term f(x) in each grid cell
  Eigen::VectorXd recharge = Eigen::VectorXd::Ones(numCells);

  // Create an instance of the FlowModel class
  auto mod = std::make_shared<FlowEquation>(recharge);

  // Define a conductivity field and evaluate the model k(x) = exp(cos(20*x))
  Eigen::VectorXd cond = (20.0*cellLocs.array()).cos().exp();

  Eigen::VectorXd h = mod->Evaluate(cond).at(0);

  std::cout << "Solution: \n" << h.transpose() << std::endl;

  /***
  ## Check Gradient of Model 
  */
  auto objective = std::make_shared<Gaussian>(numCells+1)->AsDensity();
  
  // Evaluate the gradient 
  Eigen::VectorXd objSens = Eigen::VectorXd::Ones(1);
  Eigen::VectorXd sens = objective->Gradient(0,0,h,objSens);
  Eigen::VectorXd grad = mod->Gradient(0,0,cond,sens);

  std::cout << "Gradient: \n" << grad.transpose() << std::endl;

  /***
  ## Test Hessian of Model
  */
  Eigen::VectorXd hessDir = RandomGenerator::GetUniform(numCells);

  Eigen::VectorXd hessAct = mod->ApplyHessian(0,0,0,cond,sens,hessDir);
  std::cout << "Hessian Action: \n" << hessAct.transpose() << std::endl;
  Eigen::VectorXd hessActFD = mod->ApplyHessianByFD(0,0,0,std::vector<Eigen::VectorXd>{cond},sens,hessDir);
  std::cout << " vs \n" << hessActFD.transpose() << std::endl;

  /***
  ## Construct Objective Function
  */
  WorkGraph graph;
  graph.AddNode(mod, "Model");
  graph.AddNode(objective, "Objective");
  graph.AddEdge("Model",0,"Objective",0);

  auto fullMod = graph.CreateModPiece("Objective");
    
  /***
  ## Evaluate Hessian of Objective
  */


  /***
  ## Compute Hessian Spectrum
  */

  /***
  #### Define linear operator
  */
  unsigned int outWrt = 0; // There's only one output of "fullMod", so this is the only possibility
  unsigned int inWrt1 = 0; // There's also only one input of "fullMod"
  unsigned int inWrt2 = 0;

  double scaling = -1.0;
  std::vector<Eigen::VectorXd> inputs(1);
  inputs.at(0) = cond;
  auto hessOp = std::make_shared<HessianOperator>(fullMod, inputs, outWrt, inWrt1, inWrt2, objSens, scaling);

  /***
  #### Construct eigenvalue solver
  */
  boost::property_tree::ptree opts;
  opts.put("NumEigs", numCells); // Maximum number of eigenvalues to compute
  opts.put("Verbosity", 3); // Controls how much information is printed to std::cout by the solver

  StochasticEigenSolver solver(opts);
  solver.compute(hessOp);

  Eigen::VectorXd vals = solver.eigenvalues();
  Eigen::MatrixXd vecs = solver.eigenvectors();

  std::cout << "Eigenvalues:\n" << vals.transpose() << std::endl;
  return 0;
}
