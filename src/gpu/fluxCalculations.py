'''
Created on May 21, 2013

@author: tristan
'''

import numpy as np
from pycuda.compiler import SourceModule

fluxCode = open('fluxCalculations.cu', 'r')

try:

    # Put the kernel code into a SourceModule
    fluxModule = SourceModule(fluxCode.read())
    fluxCode.close()

    # Create reference to the specific functions in the SourceModule
    FluxSolverFn = fluxModule.get_function("FluxSolver")

    # Create callable functions
    def FluxSolver(FluxesGPU, UIntPtsGPU, BottomIntPtsGPU, propSpeedsGPU, m, n, blockDims, gridDims):

        FluxSolverFn(FluxesGPU, UIntPtsGPU, BottomIntPtsGPU, propSpeedsGPU,
                     np.int32(m), np.int32(n),
                     block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

except IOError:
    print "Error opening fluxCalculations.cu"
