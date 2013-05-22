'''
Created on May 22, 2013

@author: tristan
'''

import numpy as np
from pycuda.compiler import SourceModule

sourceCode = open('sourceCalculations.cu', 'r')

try:

    # Put the kernel code into a SourceModule
    sourceModule = SourceModule(sourceCode.read())
    sourceCode.close()

    # Create reference to the specific functions in the SourceModule

    # Create callable functions

except IOError:
    print "Error opening sourceCalculations.cu"
