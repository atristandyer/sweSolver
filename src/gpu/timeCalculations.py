'''
Created on May 28, 2013

@author: tristan
'''

import os
import numpy as np
from pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray
import pycuda.cumath as cumath

timeCode = open(os.path.join(os.path.dirname(__file__), './timeCalculations.cu'), 'r')

try:

    # Put the kernel code into a SourceModule
    timeModule = SourceModule(timeCode.read())
    timeCode.close()

    # Create reference to the specific functions in the SourceModule
    UstarFn = timeModule.get_function("buildUstar")
    UnextFn = timeModule.get_function("buildUnext")

    # Create callable functions
    def buildUstar(UstarGPU, UGPU, RGPU, ShearSourceGPU, dt, m, n, blockDims, gridDims):

        UstarFn(UstarGPU, UGPU, RGPU, ShearSourceGPU,
                np.float32(dt), np.int32(m), np.int32(n),
                block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

    def buildUnext(UnextGPU, UGPU, UstarGPU, RstarGPU, ShearSourceStarGPU, dt, m, n, blockDims, gridDims):

        UnextFn(UnextGPU, UGPU, UstarGPU, RstarGPU, ShearSourceStarGPU,
                np.float32(dt), np.int32(m), np.int32(n),
                block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

except IOError:
    print "Error opening timeCalculations.cu"


# This function assumes a square cell and destructively modifies
# the PropSpeedsGPU matrix
def calculateTimestep(PropSpeedsGPU, cellDim):

    maxPropSpeed = gpuarray.max(cumath.fabs(PropSpeedsGPU)).get()
    return cellDim / (4.0 * maxPropSpeed)
