'''
Created on May 30, 2013

@author: tristan
'''

outlines = False
dirValues = {'N': 0, 'S': 1, 'E': 2, 'W': 3}
minDots = 1

def printMatrix(matrix, direction='N'):

    matrixShape = matrix.shape
    m = matrixShape[0]  # Number of rows
    n = matrixShape[1]  # Number of columns
    nDim = len(matrixShape) - 2  # Number of dimensions of data stored at each cell

    # Single value stored in each cell
    if nDim == 0:

        # Calculate the number of values to print per row
        horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 7, minDots)

        for i in range(numTopValues):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[m - i - 1][j]
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[m - i - 1][n - j - 1]
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill
        print verticalFill
        for i in reversed(range(numBottomValues)):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[i][j]
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[i][n - j - 1]
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill

    if nDim == 1:

        # Two dim array stored in each cell
        if matrixShape[2] == 2:
            horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 13, minDots)

            for i in range(numTopValues):
                print horizontalFill
                line = ''
                for j in range(numLeftValues):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f]" % (matrix[m - i - 1][j][0], matrix[m - i - 1][j][1])
                line += middleFill
                for j in reversed(range(numRightValues)):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f]" % (matrix[m - i - 1][n - j - 1][0], matrix[m - i - 1][n - j - 1][1])
                line += ('|' if outlines else ' ')
                print line
            print horizontalFill
            print verticalFill
            for i in reversed(range(numBottomValues)):
                print horizontalFill
                line = ''
                for j in range(numLeftValues):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f]" % (matrix[i][j][0], matrix[i][j][1])
                line += middleFill
                for j in reversed(range(numRightValues)):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f]" % (matrix[i][n - j - 1][0], matrix[i][n - j - 1][1])
                line += ('|' if outlines else ' ')
                print line
            print horizontalFill

        # Three dim array stored in each cell
        if matrixShape[2] == 3:
            horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 19, minDots)

            for i in range(numTopValues):
                print horizontalFill
                line = ''
                for j in range(numLeftValues):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f, %04.1f]" % (matrix[m - i - 1][j][0], matrix[m - i - 1][j][1], matrix[m - i - 1][j][2])
                line += middleFill
                for j in reversed(range(numRightValues)):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f, %04.1f]" % (matrix[m - i - 1][n - j - 1][0], matrix[m - i - 1][n - j - 1][1], matrix[m - i - 1][n - j - 1][2])
                line += ('|' if outlines else ' ')
                print line
            print horizontalFill
            print verticalFill
            for i in reversed(range(numBottomValues)):
                print horizontalFill
                line = ''
                for j in range(numLeftValues):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f, %04.1f]" % (matrix[i][j][0], matrix[i][j][1], matrix[i][j][2])
                line += middleFill
                for j in reversed(range(numRightValues)):
                    line += ('|' if outlines else ' ')
                    line += "[%04.1f, %04.1f, %04.1f]" % (matrix[i][n - j - 1][0], matrix[i][n - j - 1][1], matrix[i][n - j - 1][2])
                line += ('|' if outlines else ' ')
                print line
            print horizontalFill

        # N, S, E, W with a single value stored at each interface
        if matrixShape[2] == 4:
            dirIndex = dirValues.get(direction)
            printNSEWdirection(m, n, matrix, dirIndex, 1)

    if nDim == 2:

        # N, S, E, W matrix
        if matrixShape[2] == 4:

            dirIndex = dirValues.get(direction)
            printNSEWdirection(m, n, matrix, dirIndex, matrixShape[3])

def printNSEWdirection(m, n, matrix, dirIndex, numValuesPerDir):

    if numValuesPerDir == 1:
        horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 7, minDots)
        for i in range(numTopValues):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[m - i - 1][j][dirIndex]
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[m - i - 1][n - j - 1][dirIndex]
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill
        print verticalFill
        for i in reversed(range(numBottomValues)):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[i][j][dirIndex]
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "%06.2f" % matrix[i][n - j - 1][dirIndex]
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill
    elif numValuesPerDir == 2:
        horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 13, minDots)
        for i in range(numTopValues):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f]" % (matrix[m - i - 1][j][dirIndex][0], matrix[m - i - 1][j][dirIndex][1])
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f]" % (matrix[m - i - 1][n - j - 1][dirIndex][0], matrix[m - i - 1][n - j - 1][dirIndex][1])
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill
        print verticalFill
        for i in reversed(range(numBottomValues)):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f]" % (matrix[i][j][dirIndex][0], matrix[i][j][dirIndex][1])
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f]" % (matrix[i][n - j - 1][dirIndex][0], matrix[i][n - j - 1][dirIndex][1])
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill
    elif numValuesPerDir == 3:
        horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 19, minDots)
        for i in range(numTopValues):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f, %04.1f]" % (matrix[m - i - 1][j][dirIndex][0], matrix[m - i - 1][j][dirIndex][1], matrix[m - i - 1][j][dirIndex][2])
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f, %04.1f]" % (matrix[m - i - 1][n - j - 1][dirIndex][0], matrix[m - i - 1][n - j - 1][dirIndex][1], matrix[m - i - 1][n - j - 1][dirIndex][2])
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill
        print verticalFill
        for i in reversed(range(numBottomValues)):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f, %04.1f]" % (matrix[i][j][dirIndex][0], matrix[i][j][dirIndex][1], matrix[i][j][dirIndex][2])
            line += middleFill
            for j in reversed(range(numRightValues)):
                line += ('|' if outlines else ' ')
                line += "[%04.1f, %04.1f, %04.1f]" % (matrix[i][n - j - 1][dirIndex][0], matrix[i][n - j - 1][dirIndex][1], matrix[i][n - j - 1][dirIndex][2])
            line += ('|' if outlines else ' ')
            print line
        print horizontalFill


def getTerminalSize():
    import os
    env = os.environ
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        cr = (env.get('LINES', 25), env.get('COLUMNS', 80))

    return int(cr[1]), int(cr[0])

def createFillLines(m, n, valueWidth, minDotWidth):

    termSize = getTerminalSize()
    termWidth = termSize[0]  # The width of the terminal window
    termHeight = termSize[1]  # The height of the terminal window

    # Horizontal calculations
    extraSpace = termWidth % valueWidth
    if extraSpace < minDotWidth + 1:
        numHorizontalValues = (termWidth - extraSpace) / valueWidth - 1
        extraSpace += valueWidth
    else:
        numHorizontalValues = (termWidth - extraSpace) / valueWidth

    dotWidth = extraSpace - 1

    if numHorizontalValues % 2 == 0:
        numLeftValues = numHorizontalValues / 2
    else:
        numLeftValues = numHorizontalValues / 2 + 1
    numRightValues = numHorizontalValues / 2

    # Vertical Calculations
    numVerticalValues = termHeight / 2 - 3
    if numVerticalValues % 2 == 0:
        numTopValues = numVerticalValues / 2
    else:
        numTopValues = numVerticalValues / 2 + 1
    numBottomValues = numVerticalValues / 2

    # Build Fill Lines
    middleFill = ('+' if outlines else ' ')
    for _ in range(dotWidth - 1):
        middleFill += '.'

    horFillLine = ''
    for _ in range(numLeftValues):
        horFillLine += ('+' if outlines else ' ')
        for __ in range(valueWidth - 1):
            horFillLine += ('-' if outlines else ' ')
    horFillLine += middleFill
    for _ in range(numRightValues):
        horFillLine += ('+' if outlines else ' ')
        for __ in range(valueWidth - 1):
            horFillLine += ('-' if outlines else ' ')
    horFillLine += ('+' if outlines else ' ')

    verFillLine = ''
    for _ in range(numLeftValues):
        verFillLine += ('+' if outlines else ' ')
        for __ in range(valueWidth - 1):
            verFillLine += '.'
    verFillLine += middleFill
    for _ in range(numRightValues):
        verFillLine += ('+' if outlines else ' ')
        for __ in range(valueWidth - 1):
            verFillLine += '.'
    verFillLine += ('+' if outlines else ' ')

    middleFill = ('|' if outlines else ' ')
    for _ in range(dotWidth - 1):
        middleFill += '.'

    return horFillLine, verFillLine, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues
