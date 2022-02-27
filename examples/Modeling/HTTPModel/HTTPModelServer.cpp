#include "FlowEquation.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/HTTPModel/HTTPModPieceServer.h"

#include "MUQ/Utilities/RandomGenerator.h"

/***
## Overview

The HTTPModel interface allows coupling model and UQ codes through HTTP. A model may then
be implemented in virtually any programming language or framework, run in a container
or even on a remote machine. Likewise, the model does not make any assumptions on how the client is implemented.

This example shows how to provide a MUQ model to clients through the HTTPModel interface, specifically
the client shown in the HTTPModel Client example.
The server provides the physical model, while the client is responsible for the UQ side.

The HTTPModel interface is fully integrated in MUQ. In order to set up an HTTPModel server,
it is enough to pass a ModPiece to the serveModPiece function. All functionality provided by
the ModPiece will then be available to the client.
*/

int main(){

/***
## Set up model
First, we set up our physical model in terms of a ModPiece. We reuse the same ModPiece that is explained in the
groundwater flow equation example, and give it a simple recharge function.
*/

  unsigned int numCells = 200;
  Eigen::VectorXd recharge = Eigen::VectorXd::Ones(numCells);

  auto model = std::make_shared<FlowEquation>(recharge);

/***
## Expose model via HTTP
Providing the model as an HTTPModel via network is then just a single line. From that point on,
the server will run indefinitely, waiting for requests for model evaluations coming from clients.
Instructions for how to build a compatible client in MUQ can be found in the HTTPModel Client example.
*/

  muq::Modeling::serveModPiece(model, "0.0.0.0", 4242);

}
