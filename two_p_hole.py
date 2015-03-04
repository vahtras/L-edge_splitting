#!/usr/bin/env python
import sys
import os
import tarfile
import tempfile
sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))
from scipy.constants import alpha
import numpy
from daltools import prop
from dalmisc.sirrst import SiriusRestart
from dalmisc import twoso

SO_FACTOR = alpha**2/2
SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

def two_p_eigenvalues(targz, select_orbitals):
    tgz = tarfile.open(targz, 'r:gz')
    tgz.extractall(path=tempfile.gettempdir())

    sirius_rst = os.path.join(tempfile.gettempdir(), 'SIRIUS.RST')
    aoproper = os.path.join(tempfile.gettempdir(), 'AOPROPER')

    cmo = SiriusRestart(name=sirius_rst).cmo
    ls = get_ls(cmo, select_orbitals, aoproper)
    V = makeV(ls)
    er = get_eigen(V)
    return er
    

def get_ls(cmo, symorb, aoproper):
    orbitals = get_orbitals(cmo, symorb)
    spin_orbit_matrices = prop.read(*SO_LABELS, filename=aoproper)
    return [SO_FACTOR*orbitals.T*ls*orbitals for ls in spin_orbit_matrices]

def get_ls2(sirius_rst, symorb, ao2soint):
    rhf_density  = sirius_rst.get_rhf_density()
    #where comes the density from, i.e. nocc.
    orbitals = get_orbitals(sirius_rst.cmo, symorb)
    spin_orbit_matrices = [twoso.fock(rhf_density, c, filename=ao2soint) for c in "xyz"]
    return [SO_FACTOR*orbitals.T*ls2*orbitals for ls2 in spin_orbit_matrices]

def get_orbitals(cmo, symorb):
    indices = get_orbital_indices(cmo, symorb)
    return cmo.unblock()[:, indices]


def get_orbital_indices(cmo, symorb):
    indices = (sum(cmo.nrow[:sym-1]) + orb-1 for sym in symorb for orb in symorb[sym])
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

def main():
    
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('orbitals', help='Select orbitals as dict')
    parser.add_argument('daltargz', help='Dalton tar ball')

    args = parser.parse_args()

    exec "orbitals = %s" % args.orbitals
    print two_p_eigenvalues(args.daltargz, orbitals)

if __name__ == "__main__":
    main()
