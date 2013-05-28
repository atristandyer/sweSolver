'''
Created on May 28, 2013

@author: tristan
'''

import pycuda.autoinit
import numpy as np
import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda

def runGPUsimulation(m, n, U, Coordinates, BottomIntPts, Wind, saveInfo, runTime, boundaryConditions):

    # Tell the user we're starting the GPU simulation
    print "\n"
    print "------------------------------------------"
    print "----------Running GPU Simulation----------"
    print "------------------------------------------"
    print "\n"

    # Define the block and grid sizes
    blockDim = 16
    gridM = m / blockDim + (1 if (m % blockDim != 0) else 0)
    gridN = n / blockDim + (1 if (n % blockDim != 0) else 0)

    # Get save info
    saveOutput = saveInfo[0]
    if saveOutput:
        workingDir = saveInfo[1]
        saveTimestepInterval = saveInfo[2]

    # Initialize timestepping variables
    time = 0.0
    dt = 0.0
    nextSave = 0.0
    iterations = 0
    savedTimesteps = 0

    # Create the fort.14 file and open the fort.63g file for saving data

    # Transfer mesh data to GPU
    UGPU = sendToGPU(U)
    BottomIntPtsGPU = sendToGPU(BottomIntPts)
    WindGPU = sendToGPU(Wind)

    # Allocate other necessary arrays on GPU
    UIntPtsGPU = gpuarray.zeros((m, n, 4, 3), np.float32)
    HUVIntPtsGPU = gpuarray.zeros((m, n, 4, 3), np.float32)
    PropSpeedsGPU = gpuarray.zeros((m, n, 4), np.float32)
    FluxesGPU = gpuarray.zeros((m, n, 2, 3), np.float32)
    SlopeSourceGPU = gpuarray.zeros((m, n, 2), np.float32)
    ShearSourceGPU = gpuarray.zeros((m, n), np.float32)
    WindShearGPU = gpuarray.zeros((m, n, 2), np.float32)
    RValuesGPU = gpuarray.zeros((m, n, 3), np.float32)
    UstarGPU = gpuarray.zeros_like(UGPU)


# This is a helper function, used to send a numpy array of floats
# to the GPU
def sendToGPU(numpyArray):

    numpyArray = numpyArray.astype(np.float32)
    arrayGPU = gpuarray.to_gpu(numpyArray)
    return arrayGPU

# This function displays the current GPU memory usage
def printGPUMemUsage():
    memory = cuda.mem_get_info()
    print "Memory Usage: " + str(100 * (1 - float(memory[0]) / float(memory[1]))) + "%\t" + str(memory[0] / 1048576.0) + " MB Free"
