
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/UMBridge/UMBridgeModPiece.h"
#include "MUQ/Modeling/LinearAlgebra/HessianOperator.h"
#include "MUQ/Modeling/LinearAlgebra/StochasticEigenSolver.h"

#include <boost/property_tree/ptree.hpp>

/***
## Overview

The UM-Bridge interface allows coupling model and UQ codes through HTTP. A model may then
be implemented in virtually any programming language or framework, run in a container
or even on a remote machine. Likewise, the model does not make any assumptions on how the client is implemented.

This example shows how to connect to a running UM-Bridge server that is implemented in the UM-Bridge Server example.
The server provides the physical model, while the client is responsible for the UQ side.

The UM-Bridge interface is fully integrated in MUQ and can be used by means of the UMBridgeModPiece class.
Once such an UMBridgeModPiece is set up, it can be used like any other ModPiece. If the model supports the respective
functionality, the ModPiece then provides simple model evaluations,
gradient evaluations, applications of the Jacobian etc.

*/
int main(){

  using namespace muq::Modeling;

/***
## Connect to model server
First, we set up an UMBridgeModPiece that connects to our model server, giving the server's address
and the model name we'd like to use. This assumes that, before the client is started,
the model server from the UM-Bridge Server example is already running on your machine.
*/

  auto mod = std::make_shared<UMBridgeModPiece>("http://localhost:4242", "forward");

  std::cout << mod->inputSizes << std::endl;
  std::cout << mod->outputSizes << std::endl;

/***
## Set up dimensions
We then set up some helpers for determining dimensions and coordinates needed below.

Note that, from this point on, this example is completely identical to the FlowEquation example.
*/

  unsigned int numCells = 200;

  // Vector containing the location of each node in the "mesh"
  Eigen::VectorXd nodeLocs = Eigen::VectorXd::LinSpaced(numCells+1,0,1);

  // Vector containing the midpoint of each cell
  Eigen::VectorXd cellLocs = 0.5*(nodeLocs.head(numCells) + nodeLocs.tail(numCells));


/***
## Evaluate Model
Here we construct a vector of conductivity values and call the `FlowEquation.Evaluate` function to evaluate the model.
Recall that you should implement the `EvaluateImpl` function, but call the `Evaluate` function.
In addition to calling `EvaluateImpl`, the `Evaluate` function checks the size of the input vectors, tracks run
times, counts function evaluations, and manages a one-step cache of outputs.
*/

  // Define a conductivity field and evaluate the model k(x) = exp(cos(20*x))
  Eigen::VectorXd cond = (20.0*cellLocs.array()).cos().exp();

  Eigen::VectorXd h = mod->Evaluate(cond).at(0);

  std::cout << "Solution: \n" << h.transpose() << std::endl << std::endl;

/***
## Check Model Gradient
To check our adjoint gradient implementation, we will employ a finite difference approximation of the gradient vector.
Before doing that however, we need to define a scalar "objective function" $J(h)$ that can be composed with the flow
equation model.   In practice, this objective function is often the likelihood function or posterior density in a
Bayesian inverse problem.  For simplicity, we will just consider the log of a standard normal density:
$$
J(h) \propto -\frac{1}{2} \|h\|^2.
$$

The following cell uses the density view of MUQ's `Gaussian` class to define $J(h)$.   The `objective`
defined in this cell is just another `ModPiece` and be used in the same way as any other `ModPiece`.
*/
  auto objective = std::make_shared<Gaussian>(numCells+1)->AsDensity();

/***
To obtain the gradient of the objective function $J(h)$ with respect to the vector of cell-wise conductivities
$\mathbf{k}$, we need to apply the chain rule:
$$
\nabla_{\mathbf{k}} J = \left( \nabla_h J \right) \nabla_{\mathbf{k}} h
$$
The following cell uses the `objective` object to obtain the initial sensitivity $\nabla_h J$.  This is then
passed to the `Gradient` function of the flow model, which will use our adjoint implementation above to
compute $\nabla_{\mathbf{k}} J$.
*/
  Eigen::VectorXd objSens = Eigen::VectorXd::Ones(1);
  Eigen::VectorXd sens = objective->Gradient(0,0,h,objSens);
  Eigen::VectorXd grad = mod->Gradient(0,0,cond,sens);

  std::cout << "Gradient: \n" << grad.transpose() << std::endl << std::endl;

/***
To verify our `FlowEquation.GradientImpl` function, we can call the built in `ModPiece.GradientByFD` function
to construct a finite difference approximation of the gradient.  If all is well, the finite difference and adjoint
gradients will be close.
*/
  Eigen::VectorXd gradFD = mod->GradientByFD(0,0,std::vector<Eigen::VectorXd>{cond},sens);
  std::cout << "Finite Difference Gradient:\n" << gradFD.transpose() << std::endl << std::endl;

/***
## Test Jacobian of Model
Here we randomly choose a vector $v$ (`jacDir`) and compute the action of the Jacobian $Jv$ using both our adjoint method and MUQ's built-in finite difference implementation.
*/
  Eigen::VectorXd jacDir = Eigen::VectorXd::Ones(numCells);//RandomGenerator::GetUniform(numCells);

  Eigen::VectorXd jacAct = mod->ApplyJacobian(0,0, cond, jacDir);
  std::cout << "Jacobian Action: \n" << jacAct.transpose() << std::endl << std::endl;

  Eigen::VectorXd jacActFD = mod->ApplyJacobianByFD(0,0, std::vector<Eigen::VectorXd>{cond}, jacDir);
  std::cout << "Finite Difference Jacobian Action \n" << jacActFD.transpose() << std::endl << std::endl;

/***
## Test Hessian of Model
We now take a similar approach to verifying our Hessian action implementation.  Here we randomly choose a vector $v$ (`hessDir`)
and compute $Hv$ using both our adjoint method and MUQ's built-in finite difference implementation.
*/
  Eigen::VectorXd hessDir = Eigen::VectorXd::Ones(numCells);//RandomGenerator::GetUniform(numCells);

  Eigen::VectorXd hessAct = mod->ApplyHessian(0,0,0,cond,sens,hessDir);
  std::cout << "Hessian Action: \n" << hessAct.transpose() << std::endl << std::endl;

  Eigen::VectorXd hessActFD = mod->ApplyHessianByFD(0,0,0,std::vector<Eigen::VectorXd>{cond},sens,hessDir);
  std::cout << "Finite Difference Hessian Action \n" << hessActFD.transpose() << std::endl << std::endl;

/***
## Test Hessian of Objective
In the tests above, we manually evaluate the `objective` and `mod` components separately.   They can also be combined
in a MUQ `WorkGraph`, which is more convenient when a large number of components are used or Hessian information needs
to be propagated through multiple different components.

The following code creates a `WorkGraph` that maps the output of the flow model to the input of the objective function.
It then creates a new `ModPiece` called `fullMod` that evaluates the composition $J(h(k))$.
*/
  WorkGraph graph;
  graph.AddNode(mod, "Model");
  graph.AddNode(objective, "Objective");
  graph.AddEdge("Model",0,"Objective",0);

  auto fullMod = graph.CreateModPiece("Objective");

/***
As before, we can apply the Hessian of the full model to the randomly generated `hessDir` vector and compare the results
with finite differences.   Notice that the results shown here are slightly different than the Hessian actions computed above.
Above, we manually fixed the sensitivity $s$ independently of $h$ and did not account for the relationship between the
conductivity $k$ on the sensitivity $s$.   The `WorkGraph` however, captures all of those dependencies.
*/

  hessAct = fullMod->ApplyHessian(0,0,0,cond,objSens,hessDir);
  hessActFD = fullMod->ApplyHessianByFD(0,0,0,std::vector<Eigen::VectorXd>{cond},objSens,hessDir);

  std::cout << "Hessian Action: \n" << hessAct.transpose() << std::endl << std::endl;

  std::cout << "Finite Difference Hessian Action \n" << hessActFD.transpose() << std::endl << std::endl;

  return 0;
}
