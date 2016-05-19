#!/usr/bin/env python
import sys
import os
import tarfile
import tempfile
sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))
from scipy.constants import alpha
import numpy
from daltools import prop
from daltools.sirrst import SiriusRestart
from two import twoso

SO_FACTOR = alpha**2/2
SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

def two_p_eigenvalues(targz, select_orbitals, two_electron=False, all_electron=False):
    """
    Calculate eigenvalues for projected spin-orbit Hamiltonian
    """

    sirius_rst, aoproper =  \
        unpack_dalfiles(targz, getfiles=['SIRIUS.RST', 'AOPROPER'])

    restart_file = SiriusRestart(name=sirius_rst)
    cmo = restart_file.cmo

    if two_electron:
        ao2soint, = unpack_dalfiles(targz, getfiles=['AO2SOINT'])
        ls = get_ls2(restart_file, select_orbitals, ao2soint)
    elif all_electron:
        ao2soint, = unpack_dalfiles(targz, getfiles=['AO2SOINT'])
        ls1 = get_ls1(cmo, select_orbitals, aoproper)
        ls2 = get_ls2(restart_file, select_orbitals, ao2soint)
        ls = [m1 + m2 for m1, m2 in zip(ls1, ls2)]
    else:
        ls = get_ls1(cmo, select_orbitals, aoproper)
    V = makeV(ls)
    er = get_eigen(V)
    return er

def unpack_dalfiles(targz, getfiles=[]):
    tgz = tarfile.open(targz, 'r:gz')
    tgz.extractall(path=tempfile.gettempdir())
    return  (os.path.join(tempfile.gettempdir(), f) for f in getfiles)
    

def get_ls1(cmo, symorb, aoproper):
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
    parser.add_argument('--two-electron', action='store_true', help='Only two-electron spin-orbit')
    parser.add_argument('--all-electron', action='store_true', help='Full Breit-Pauli spin-orbit')

    args = parser.parse_args()

    exec "orbitals = %s" % args.orbitals
    print two_p_eigenvalues(args.daltargz, orbitals, two_electron=args.two_electron, all_electron=args.all_electron)

if __name__ == "__main__":
    main()
