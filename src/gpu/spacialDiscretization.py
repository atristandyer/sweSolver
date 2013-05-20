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
    spacialModule = SourceModule(spacialCode.read())
except IOError:
    print "Error opening spacialDiscretization.cu"
