# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import sys
import numpy

sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))

# <codecell>

SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')
eV = 1/27.211396132

# <codecell>

def two_p_splitting(cmo, occa, occb, ps, aoproper="AOPROPER", ao2soint="AO2SOINT", one_electron=False):
    
    SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')
    
    from scipy.constants import alpha
    prefactor = alpha**2/2
    
    occupied_alpha = cmo.get_columns(occa)
    occupied_beta = cmo.get_columns(occb)
    
    Da = (occupied_alpha*occupied_alpha.T()).unblock()
    Db = (occupied_beta*occupied_alpha.T()).unblock()
    
        
    ls1 = prop.read(*SO_LABELS, filename=aoproper)
    
    if one_electron:
        ls = ls1
    else:
        ls2 = [0.5*(fa - fb) for fa, fb in [twoso.fockab(Da, Db, c, filename=ao2soint) for c in "xyz"]]  
        ls = [one + two for one, two in zip(ls1, ls2)]
    
    p = cmo.unblock()[:, ps]
    plsp = [prefactor*p.T*ls_component*p for ls_component in ls]
    
    V = makeV(plsp)
    
    eigvals = numpy.linalg.eigvals(V)
    er = numpy.sort(eigvals.real)
    
    print er
    split = abs(.5*er[:2].sum() - .25*er[2:].sum())
    
    return split/eV

def makeV(ls):
    V = numpy.zeros((6, 6), dtype='complex', order='F')
    V[:3, :3] = ls[2]
    V[3:, 3:] = -ls[2]
    V[:3, 3:] = ls[0] + 1j*ls[1]
    V[3:, :3] = ls[0] - 1j*ls[1]
    V *= 1j/2
    return V


# <codecell>

from dalmisc import sirrst, twoso
from daltools import prop

# <codecell>

#system specific
mol = ("Cystein", "cys")
mol = ("S", "S")
aoproper = "%s/hf_%s.AOPROPER" % mol
ao2soint = "%s/hf_%s.AO2SOINT" % mol
sirius_rst = "%s/hf_%s.SIRIUS.RST" % mol

nocc_alpha = [32] #[3, 2, 2, 0, 1, 0, 0, 0]
nocc_beta  = [32] #[3, 2, 2, 0, 1, 0, 0, 0]

# <codecell>

!ls Cystein

# <codecell>


cmo = sirrst.SiriusRestart(name=sirius_rst).cmo

# <codecell>

x, y, z = 9, 15, 23
#x, y, z = 9, 10, 11 #9, 15, 23

# <codecell>

print nocc_alpha, nocc_beta

# <codecell>

print cmo.nrow

# <codecell>

two_p_splitting(cmo, nocc_alpha, nocc_beta, (x, y, z), aoproper=aoproper, ao2soint=ao2soint)

# <codecell>

cmorel = sirrst.SiriusRestart(name="S/hfrel_S.SIRIUS.RST").cmo

# <codecell>

two_p_splitting(cmorel, nocc_alpha, nocc_beta, (x, y, z), aoproper=aoproper, ao2soint=ao2soint, one_electron=True)

# <codecell>


