import sys, os
sys.path.append(os.environ['MUQ_DIR']+'/../../../'+'python')

import muq.Modeling as mm
import muq.Utilities as mu
import muq.Approximation as ma

import numpy as np

logStiffMean = [8.0, 8.0, 8.0]
logStiffStd = [1.5, 1.5, 1.5]

fub = -0.0
flb = -0.1
loadMean = [-0.05]
loadScale = [0.05]

stiffTransform = mm.AffineOperator(mm.DiagonalOperator(logStiffStd), logStiffMean)
loadTransform = mm.AffineOperator(mm.DiagonalOperator([0.5*(fub-flb)]), [0.5*(fub+flb)])

print(stiffTransform)
