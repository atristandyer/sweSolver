'''
Created on May 20, 2013

@author: tristan
'''

import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

spacialCode = open('spacialDiscretization.cu', 'r')

try:

    # Put the kernel code into a SourceModule
    spacialModule = SourceModule(spacialCode.read())
    spacialCode.close()

    # Create reference to the specific functions in the SourceModule
    ReconstructFreeSurfaceFn = spacialModule.get_function("ReconstructFreeSurface")

    # Create callable functions
    def ReconstructFreeSurface(UGPU, BottomIntPtsGPU, UIntPtsGPU, huvIntPtsGPU, m, n, dx, dy, blockDims, gridDims):

        ReconstructFreeSurfaceFn(UGPU, BottomIntPtsGPU, UIntPtsGPU, huvIntPtsGPU,
                                 np.int32(m), np.int32(n), np.float32(dx), np.float32(dy),
                                 block=(blockDims[0], blockDims[1], 1), grid=(gridDims[0], gridDims[1]))

except IOError:
    print "Error opening spacialDiscretization.cu"
