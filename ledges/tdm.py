"""Transitions for in general non-orthogonal orbitals"""
from vb.nod import Nod, ao_transition_matrix
import util

def get_transition_spin_densities(*nods):
    """Calculate matrix of spin_density matrices"""
    nbas = Nod.S.shape[0]
    dims = [nbas] * len(nods)
    dens = util.subblocked.matrix(dims, dims)
    for i, nod_i in enumerate(nods):
        for j, nod_j in enumerate(nods):
            d_a, d_b = ao_transition_matrix(nod_i, nod_j)
            dens.subblock[i][j] = (d_a - d_b)/2
    return dens
