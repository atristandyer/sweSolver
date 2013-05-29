'''
Created on May 29, 2013

@author: tristan
'''

import numpy as np
import math

# This functions builds a calm, flat initial water surface with constant elevation
def buildFlatWaterSurface(m, n, elevation):

    meshU = np.zeros((m, n, 3))

    for i in range(m):
        for j in range(n):
            meshU[i][j][0] = elevation

    return meshU


# This function builds a flat initial water surface with a water pyramid located at the given
# x, y coordintes. It uses the iterative pyramid builder along with the flat water surface builder
def buildPyramidWaterSurface(m, n, elevation, x, y, pyrElevation, dz):

    radius = int(math.floor((pyrElevation - elevation) / dz))
    meshU = buildPyramid(buildFlatWaterSurface(m, n, elevation), m, n, x, y, radius, dz)
    return meshU

# Recursively builds a water pyramid centered at centerX, centerY by
# adding the value of height to the surface at each iteration
def buildPyramid(meshU, m, n, centerX, centerY, radius, height):

    if (radius == 0):
        meshU[centerY][centerX][0] += height
        return meshU
    else:
        for i in range(centerY - radius, centerY + radius):
            for j in range(centerX - radius, centerX + radius):
                if (i >= 0 and i < m and j >= 0 and j < n):
                    meshU[i][j][0] += height
        return buildPyramid(meshU, m, n, centerX, centerY, radius - 1, height)

# This function is used to create the wind source matrix. A constant value could be used
# instead of building an entire matrix with the same values, but this will eventually
# allow for dynamic wind conditions.
def buildWindSource(m, n, speed, direction):

    meshWind = np.zeros((m, n, 2))
    for i in range(m):
        for j in range(n):
            meshWind[i][j][0] = speed
            meshWind[i][j][1] = math.radians(direction)

    return meshWind

def buildLinearWindSource(m, n, speed, direction):

    meshWind = np.zeros((m, n, 2))
    for i in range(m):
        for j in range(n):
            meshWind[i][j][0] = (2.0 / 3.0) * speed - (float(j) / float(n)) * speed
            meshWind[i][j][1] = math.radians(direction)

    return meshWind


# This function is used to reset the water elevations for all dry areas
# before beginning a simulation
def validateInitialConditions(m, n, meshBottomIntPts, meshU):

    for i in range(m):
        for j in range(n):
            cellCenter = (meshBottomIntPts[i][j + 1][1] + meshBottomIntPts[i][j][1]) / 2.0
            if meshU[i][j][0] < cellCenter:
                meshU[i][j][0] = cellCenter
    return meshU
