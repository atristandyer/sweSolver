'''
Created on May 29, 2013

@author: tristan
'''


################ Dependency Tests #################
try:
    import numpy as np
except ImportError:
    print "Please install Numpy"
    exit(0)

try:
    import pycuda.autoinit
    import pycuda.driver as cuda
    gpuValid = True
except ImportError:
    gpuValid = False
    print "Unable to load PyCUDA. Please install it for GPU functionality"

try:
    from mesh.meshBuilder import buildSlopingDomain
    from gpu.gpuSimulation import runGPUsimulation
except ImportError:
    print "Error loading sweSolver libraries, please reinstall"
    exit(0)



######################################################
##### Testing
######################################################
Coordinates, BottomIntPts = buildSlopingDomain(1.0, 128, 128, 0.0, 0.0, 1.0, 0)
print Coordinates.shape
print BottomIntPts.shape
