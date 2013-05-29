'''
Created on May 20, 2013

@author: tristan
'''

import os
import numpy as np
from pycuda.compiler import SourceModule
from gpu.gpuSimulation import useCachedKernels

spacialCode = open(os.path.join(os.path.dirname(__file__), 'spacialDiscretization.cu'), 'r')

try:

    # Put the kernel code into a SourceModule
    if useCachedKernels:
        spacialModule = SourceModule(spacialCode.read())
    else:
        spacialModule = SourceModule(spacialCode.read(), cache_dir=False)
    spacialCode.close()

    # Create reference to the specific functions in the SourceModule
    ReconstructFreeSurfaceFn = spacialModule.get_function("ReconstructFreeSurface")
    CalculatePropSpeedsFn = spacialModule.get_function("CalculatePropSpeeds")

    # Create callable functions
    def ReconstructFreeSurface(UGPU, BottomIntPtsGPU, UIntPtsGPU, huvIntPtsGPU, m, n, dx, dy, blockDims, gridDims):

        ReconstructFreeSurfaceFn(UGPU, BottomIntPtsGPU, UIntPtsGPU, huvIntPtsGPU,
                                 np.int32(m), np.int32(n), np.float32(dx), np.float32(dy),
                                 block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

    def CalculatePropSpeeds(UIntPtsGPU, huvIntPtsGPU, propSpeedsGPU, m, n, blockDims, gridDims):

        CalculatePropSpeedsFn(UIntPtsGPU, huvIntPtsGPU, propSpeedsGPU,
                              np.int32(m), np.int32(n),
                              block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

except IOError:
    print "Error opening spacialDiscretization.cu"
