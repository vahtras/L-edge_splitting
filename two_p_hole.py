#!/usr/bin/env python
import sys
import os
import tarfile
import tempfile
sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))
from scipy.constants import alpha
import numpy
from daltools import prop

SO_FACTOR = alpha**2/2
SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

def two_p_eigenvalues(targz, select_orbitals):
    tgz = tarfile.open(targz, 'r:gz')
    tgz.extractall(tempfile.gettempdir(), ['SIRIUS.RST', 'AOPROPER'])

def get_ls(cmo, symorb, aoproper):
    orbitals = get_orbitals(cmo, symorb)
    spin_orbit_matrices = prop.read(*SO_LABELS, filename=aoproper)
    return [SO_FACTOR*orbitals.T*ls*orbitals for ls in spin_orbit_matrices]

def get_orbitals(cmo, symorb):
    indices = get_orbital_indices(cmo, symorb)
    return cmo.unblock()[:, (9, 15, 23)]


def get_orbital_indices(cmo, symorb):
    indices = (sum(cmo.nrow[:sym-1]) for sym in symorb)
    return tuple(indices)

    
def makeV(ls):
    V = numpy.zeros((6, 6), dtype='complex', order='F')
    V[:3, :3] = ls[2]
    V[3:, 3:] = -ls[2]
    V[:3, 3:] = ls[0] + 1j*ls[1]
    V[3:, :3] = ls[0] - 1j*ls[1]
    V *= 1j/2
    return V

def get_eigen(V):
    eigenvalues = numpy.linalg.eigvals(V)
    return numpy.sort(eigenvalues.real)
