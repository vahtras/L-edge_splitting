#!/usr/bin/env python
import sys
import os
sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))
from daltools import prop
from scipy.constants import alpha
SO_FACTOR = alpha**2/2

SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

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

    
