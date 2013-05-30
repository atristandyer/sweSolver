'''
Created on May 30, 2013

@author: tristan
'''

outlines = False

def printMatrix(matrix, direction='N'):

    matrixShape = matrix.shape
    m = matrixShape[0]  # Number of rows
    n = matrixShape[1]  # Number of columns
    nDim = len(matrixShape) - 2  # Number of dimensions of data stored at each cell

    # Single value stored in each cell
    if nDim == 0:

        # Calculate the number of values to print per row
        horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 7, 2)

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

    # Data stored in each cell is an array
    if 0 == 1:

        # Calculate the number of values to print per row
        horizontalFill, verticalFill, middleFill, numLeftValues, numRightValues, numTopValues, numBottomValues = createFillLines(m, n, 7, 2)

        for i in range(numTopValues):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += '|      '
            line += middleFill
            for j in range(numRightValues):
                line += '|      '
            line += '|'
            print line
        print horizontalFill
        print verticalFill
        for i in range(numBottomValues):
            print horizontalFill
            line = ''
            for j in range(numLeftValues):
                line += '|      '
            line += middleFill
            for j in range(numRightValues):
                line += '|      '
            line += '|'
            print line
        print horizontalFill


        # Single value stored in each cell
        if matrixShape[2] == 1:
            for i in reversed(range(n - numTopValues, n)):
                line = '|'
                for j in range(0, numLeftValues):
                    line += "%05.2f" % matrix[i][j] + ' | '
                print line


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
