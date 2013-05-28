'''
Created on May 22, 2013

@author: tristan
'''

import numpy as np
from pycuda.compiler import SourceModule

sourceCode = open('sourceCalculations.cu', 'r')

try:

    # Put the kernel code into a SourceModule
    sourceModule = SourceModule(sourceCode.read())
    sourceCode.close()

    # Create reference to the specific functions in the SourceModule
    BedSlopeFn = sourceModule.get_function("bedSlopeSourceSolver")
    BedShearFn = sourceModule.get_function("bedShearSourceSolver")

    # Create callable functions
    def BedSlopeSourceSolver(SlopeSourceGPU, UGPU, BottomIntPtsGPU, m, n, dx, dy, blockDims, gridDims):

        BedSlopeFn(SlopeSourceGPU, UGPU, BottomIntPtsGPU,
                   np.int32(m), np.int32(n), np.float32(dx), np.float32(dy),
                   block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

    def BedShearSourceSolver(ShearSourceGPU, UGPU, BottomIntPtsGPU, m, n, dx, dy, blockDims, gridDims):

        BedShearFn(ShearSourceGPU, UGPU, BottomIntPtsGPU,
                   np.int32(m), np.int32(n), np.float32(dx), np.float32(dy),
                   block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

except IOError:
    print "Error opening sourceCalculations.cu"
