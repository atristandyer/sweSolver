'''
Created on May 28, 2013

@author: tristan
'''

import numpy as np
from pycuda.compiler import SourceModule

boundaryCode = open('boundaryConditions.cu', 'r')

try:

    # Put the kernel code into a SourceModule
    boundaryModule = SourceModule(boundaryCode.read())
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
