'''
Created on May 28, 2013

@author: tristan
'''

from sys import stdout
import pycuda.autoinit
import numpy as np
import time as timer
import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda

# Import the kernels used for performing calculations on the GPU
from spacialDiscretization import ReconstructFreeSurface, CalculatePropSpeeds
from fluxCalculations import FluxSolver, BuildRValues
from sourceCalculations import BedSlopeSourceSolver, BedShearSourceSolver
from timeCalculations import calculateTimestep, buildUstar, buildUnext
from boundaryConditions import ApplyWallBoundaries

def runGPUsimulation(m, n, U, Coordinates, BottomIntPts, Wind, saveInfo, runTime, boundaryConditions):

    # Show GPU memory usage before starting the simulation
    printGPUMemUsage()

    # Tell the user we're starting the GPU simulation
    print "\n"
    print "------------------------------------------"
    print "----------Running GPU Simulation----------"
    print "------------------------------------------"
    print "\n"

    # Calculate the cell size from the bottom left cell
    dx = Coordinates[0][1][0] - Coordinates[0][0][0]
    dy = Coordinates[1][0][1] - Coordinates[0][0][1]

    # Define the block and grid sizes
    blockDim = 16
    gridM = m / blockDim + (1 if (m % blockDim != 0) else 0)
    gridN = n / blockDim + (1 if (n % blockDim != 0) else 0)
    blockDims = [blockDim, blockDim]
    gridDims = [gridN, gridM]

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
    WindSourceGPU = gpuarray.zeros((m, n, 2), np.float32)
    RGPU = gpuarray.zeros((m, n, 3), np.float32)
    UstarGPU = gpuarray.zeros_like(UGPU)

    # Print the GPU's memory usage
    printGPUMemUsage()

    # Begin timestepping
    sTime = timer.time()
    while time < runTime:

        # Save output

        ###################################################
        ##### Begin first order accurate calculations #####
        ###################################################

        # Reconstruct the free surface
        ReconstructFreeSurface(UGPU, BottomIntPtsGPU, UIntPtsGPU, HUVIntPtsGPU, m, n, dx, dy, blockDims, gridDims)

        # Calculate propagation speeds
        CalculatePropSpeeds(UIntPtsGPU, HUVIntPtsGPU, PropSpeedsGPU, m, n, blockDims, gridDims)

        # Calculate fluxes
        FluxSolver(FluxesGPU, UIntPtsGPU, BottomIntPtsGPU, PropSpeedsGPU, m, n, blockDims, gridDims)

        # Calculate source terms
        BedSlopeSourceSolver(SlopeSourceGPU, UGPU, BottomIntPtsGPU, m, n, dx, dy, blockDims, gridDims)
        BedShearSourceSolver(ShearSourceGPU, UGPU, BottomIntPtsGPU, m, n, dx, dy, blockDims, gridDims)

        # Build R
        BuildRValues(RGPU, FluxesGPU, SlopeSourceGPU, WindSourceGPU, m, n, blockDims, gridDims)

        # Calculate timestep
        dt = calculateTimestep(PropSpeedsGPU, dx)
        if saveOutput and time + dt > nextSave:
            dt = nextSave - time

        # Build U*
        buildUstar(UstarGPU, UGPU, RGPU, ShearSourceGPU, dt, m, n, blockDims, gridDims)

        # Apply boundary conditions
        ApplyWallBoundaries(UstarGPU, m, n, blockDims, gridDims)


        ####################################################
        ##### Begin second order accurate calculations #####
        ####################################################

        # Reconstruct the free surface
        ReconstructFreeSurface(UstarGPU, BottomIntPtsGPU, UIntPtsGPU, HUVIntPtsGPU, m, n, dx, dy, blockDims, gridDims)

        # Calculate propagation speeds
        CalculatePropSpeeds(UIntPtsGPU, HUVIntPtsGPU, PropSpeedsGPU, m, n, blockDims, gridDims)

        # Calculate fluxes
        FluxSolver(FluxesGPU, UIntPtsGPU, BottomIntPtsGPU, PropSpeedsGPU, m, n, blockDims, gridDims)

        # Calculate source terms
        BedSlopeSourceSolver(SlopeSourceGPU, UstarGPU, BottomIntPtsGPU, m, n, dx, dy, blockDims, gridDims)
        BedShearSourceSolver(ShearSourceGPU, UstarGPU, BottomIntPtsGPU, m, n, dx, dy, blockDims, gridDims)

        # Build R*
        BuildRValues(RGPU, FluxesGPU, SlopeSourceGPU, WindSourceGPU, m, n, blockDims, gridDims)

        # Build Unext
        buildUnext(UGPU, UGPU, UstarGPU, RGPU, ShearSourceGPU, dt, m, n, blockDims, gridDims)

        # Apply boundary conditions
        ApplyWallBoundaries(UGPU, m, n, blockDims, gridDims)


        ################################
        ##### Advance the timestep #####
        ################################

        time += dt
        iterations += 1

        # Print some output to the user every 100 iterations
        if (iterations % 100 == 0):
            stdout.write("\rIteration: %i\tTotal time simulated: %.4f seconds\tTimestep: %.4f" % (iterations, time, dt))
            stdout.flush()

    fTime = timer.time()

    # Save output

    # Print run information

    return fTime - sTime


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
