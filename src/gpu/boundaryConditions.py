'''
Created on May 28, 2013

@author: tristan
'''

import os
import numpy as np
from pycuda.compiler import SourceModule
from gpu.gpuSimulation import useCachedKernels

boundaryCode = open(os.path.join(os.path.dirname(__file__), './boundaryConditions.cu'), 'r')

try:

    # Put the kernel code into a SourceModule
    if useCachedKernels:
        boundaryModule = SourceModule(boundaryCode.read())
    else:
        boundaryModule = SourceModule(boundaryCode.read(), cache_dir=False)
    boundaryCode.close()

    # Create reference to the specific functions in the SourceModule
    WallFn = boundaryModule.get_function("applyWallBoundaries")

    # Create callable functions
    def ApplyWallBoundaries(UGPU, m, n, blockDims, gridDims):

        WallFn(UGPU,
               np.int32(m), np.int32(n),
               block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

except IOError:
    print "Error opening boundaryConditions.cu"
