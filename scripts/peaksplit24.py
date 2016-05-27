#!/usr/bin/env python
import sys
import numpy
from scipy.constants import physical_constants as pc
eV = pc['Hartree energy in eV'][0]

eigenvalues = numpy.loadtxt(sys.argv[1])
split = sum(eigenvalues[2:])/4 - sum(eigenvalues[:2])/2
print "Mean peak split", split
print "Mean peak split (eV)", split*eV
