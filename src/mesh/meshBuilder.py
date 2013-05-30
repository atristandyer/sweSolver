'''
Created on May 29, 2013

@author: tristan
'''


import numpy as np


# This function is used to build a sloping domain
#     cellSize - the horizontal dimension of the square cell
#     m - the number of rows
#     n - the number of columns
#     bL - the elevation of the bottom left point
#     tL - the elevation of the top left point
#     tR - the elevation of the top right point
#     boundaryType - the boundary type to be used for the entire domain
#         0 - Wall boundaries
def buildSlopingDomain(cellSize, m, n, bL, tL, tR, boundaryType):

    slopeX = (tR - tL) / (n * cellSize)
    slopeY = (tL - bL) / (m * cellSize)

    # Extra row and column needed to full define an mxn grid
    meshCoordinates = np.zeros((m + 1, n + 1, 3))
    meshBottomIntPts = np.zeros((m + 1, n + 1, 2))

    # Create all initial values
    for i in range(m + 1):
        rowHeadElevation = bL + i * slopeY
        for j in range(n + 1):
            # x-coordinate
            meshCoordinates[i][j][0] = j * cellSize
            # y-coordinate
            meshCoordinates[i][j][1] = i * cellSize
            # z-coordinate (measured from z = 0.0)
            meshCoordinates[i][j][2] = rowHeadElevation + j * slopeX

    # Change boundary elevations
    adjustBoundaries(m, n, meshCoordinates, boundaryType)

    for i in range(m + 1):
        for j in range(n + 1):
            if j < n:
                # Bottom integration point elevation
                meshBottomIntPts[i][j][0] = (meshCoordinates[i][j + 1][2] + meshCoordinates[i][j][2]) / 2.0
            if i < m:
                # Left integration point elevation
                meshBottomIntPts[i][j][1] = (meshCoordinates[i + 1][j][2] + meshCoordinates[i][j][2]) / 2.0

    return meshCoordinates, meshBottomIntPts


def adjustBoundaries(m, n, meshCoordinates, boundaryType):

    # Wall Boundaries
    if boundaryType == 0:

        # Reflect bottom row
        for i in range(2):
            for j in range(n + 1):
                meshCoordinates[i][j][2] = meshCoordinates[4 - i][j][2]

        # Reflect top row
        for i in range(m - 1, m + 1):
            for j in range(n + 1):
                meshCoordinates[i][j][2] = meshCoordinates[2 * m - 6 - i][j][2]

        # Reflect left column
        for i in range(n):
            for j in range(2):
                meshCoordinates[i][j][2] = meshCoordinates[i][4 - j][2]

        # Reflect right column:
        for i in range(n):
            for j in range(n - 1, n + 1):
                meshCoordinates[i][j][2] = meshCoordinates[i][2 * n - 6 - j][2]

    return meshCoordinates

