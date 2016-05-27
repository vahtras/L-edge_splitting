#!/usr/bin/env python
import sys
import numpy
from scipy.constants import physical_constants as pc
eV = pc['Hartree energy in eV'][0]

def mean(n1, n2, filename):
    eigenvalues = numpy.loadtxt(filename)
    assert n1 + n2 == eigenvalues.shape[0]
    split = abs(sum(eigenvalues[:n1])/n1 - sum(eigenvalues[n1:])/(n2))
    return split
