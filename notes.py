# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os
import numpy

sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))

# <markdowncell>

# **One-electron**
# 
# We need to evaluate $\langle p_{i\beta}| l^z s^z | p_{j\beta}\rangle 
# =  \langle a_{i\beta}^\dagger l^z_{pq} s^z_{pq} a_{j\beta}\rangle
# =  -\langle  {l_z}_{ji} (-1/2) \rangle 
# =  \frac 1 2 l^z_{ji}$

# <markdowncell>

# Since spin-coupling is now same for all components

# <markdowncell>

# $ V_{ij} = \frac 1 2 \begin{pmatrix}  
# l^z_{ji} & (l^x + i l^y)_{ji} \\
# (l^x + i l^y)_{ji} & -l^z_{ji}
# \end{pmatrix}$

# <markdowncell>

# **Sulphur atom**

# <markdowncell>

# Run Dalton

# <codecell>

#!cd S && ~/dev/dalton/git/build/master/dalton -get "AOONEINT AOPROPER SIRIUS.RST AO2SOINT" hf S

# <markdowncell>

# mostly we will use integrals from this set. Define filenames

# <codecell>

aooneint = "S/hf_S.AOONEINT"
aoproper = "S/hf_S.AOPROPER"
ao2soint = "S/hf_S.AO2SOINT"

# <markdowncell>

# The default input gives singlet (3px^2 3py^2 - HF occ 32201000). 
# This is triplet Sulphur (3px 3py 3pz^2) - HF occ 31102000/01100000

# <codecell>

#!cd S && ~/dev/dalton/git/build/master/dalton -get "SIRIUS.RST" hf3 S

# <markdowncell>

# Another case is relaxation, consider the ionized form(3px^2 3py^2 2pz^1)-HF occ 32200000/32201000 (for x, y we need to do inner shell hole)

# <codecell>

#!cd S && ~/dev/dalton/git/build/master/dalton -get "SIRIUS.RST" hfrel S

# <markdowncell>

# Now the integral is in S/hf_S.AOPROPER

# <codecell>

from daltools import prop
from dalmisc import sirrst

# <codecell>

angmom = prop.read('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB', filename=aoproper)

# <markdowncell>

# Get the p-orbitals, read from SIRRST

# <codecell>

rst = sirrst.SiriusRestart(name="S/hf_S.SIRIUS.RST")
cmo = rst.cmo

# <codecell>


rst3 = sirrst.SiriusRestart(name="S/hf3_S.SIRIUS.RST")
cmo3 = rst3.cmo

# <codecell>

rstrel = sirrst.SiriusRestart(name="S/hfrel_S.SIRIUS.RST")
cmorel = rstrel.cmo

# <codecell>

def get_ps(mo):
    _px, _py, _pz = (1, 2, 4)
    return mo[_px][:, 0], mo[_py][:, 0], mo[_pz][:, 0]

px, py, pz = get_ps(cmo)
px3, py3, pz3 = get_ps(cmo3)
pxrel, pyrel, pzrel = get_ps(cmorel)

#print px, py, pz
#print px3, py3, pz3
#print pxrel, pyrel, pzrel

# <codecell>

cmo0 = cmo.unblock()
cmo03 = cmo3.unblock()
cmo0rel = cmorel.unblock()

# <codecell>

x, y, z = 9, 15, 23
p = cmo0[:, (x, y, z)]
p3 = cmo03[:, (x, y, z)]
prel = cmo0rel[:, (x, y, z)]

# <codecell>

#print p, p3, prel

# <codecell>

from scipy.constants import alpha
prefactor = alpha**2/2
ls = [prefactor*p.T*M*p for M in angmom]
ls3 = [prefactor*p3.T*M*p3 for M in angmom]
lsrel = [prefactor*prel.T*M*prel for M in angmom]

# <markdowncell>

# Form complex matrices

# <codecell>

import numpy as np

# <codecell>

V = np.zeros((6, 6), dtype='complex', order='F')
V3 = np.zeros((6, 6), dtype='complex', order='F')
Vrel = np.zeros((6, 6), dtype='complex', order='F')

# <codecell>

def makeV(ls):
    V = np.zeros((6, 6), dtype='complex', order='F')
    V[:3, :3] = ls[2]
    V[3:, 3:] = -ls[2]
    V[:3, 3:] = ls[0] + 1j*ls[1]
    V[3:, :3] = ls[0] - 1j*ls[1]
    V *= 1j/2
    return V


V = makeV(ls)
V3 = makeV(ls3)
Vrel = makeV(lsrel)

# <codecell>

eigvals =  numpy.linalg.eigvals(V)
eigvals3 =  numpy.linalg.eigvals(V3)
eigvalsrel =  numpy.linalg.eigvals(Vrel)

# <codecell>

print eigvals.real
print eigvals3.real
print eigvalsrel.real

# <codecell>

eV = 27.211396132

# <codecell>

split = abs(eigvals[0] - .5*(eigvals[-1]+eigvals[-2]))
print split*eV, 'eV'

# <codecell>

split = abs(eigvals3[0] - .5*(eigvals3[-1]+eigvals3[-2]))
print split*eV, 'eV'

# <codecell>

split = abs(eigvalsrel[0] - .5*(eigvalsrel[-1]+eigvalsrel[2]))
print split*eV, 'eV'

# <markdowncell>

# **Two-electron**

# <codecell>

from dalmisc.twoso import fockab

# <markdowncell>

# The fock routine takes alpha- and beta densities

# <codecell>

#print cmo

# <codecell>

nocc_alpha = [3, 2, 2, 0, 1, 0, 0, 0]
nocc_beta  = [3, 2, 2, 0, 1, 0, 0, 0]
cmoa = [ c[:, :na] for c, na in zip(cmo, nocc_alpha) ]

# <codecell>

#for ca in cmoa:
#    print ca

# <codecell>

from daltools import dens #is there a blocked form?

# <markdowncell>

# Is there a blocked range selector for occupied orbitals?

# <codecell>

cmo_occ_alpha = cmo.get_columns(nocc_alpha)
cmo_occ_beta = cmo.get_columns(nocc_beta)

# <codecell>

Da = (cmo_occ_alpha*cmo_occ_alpha.T()).unblock()
Db = (cmo_occ_alpha*cmo_occ_alpha.T()).unblock()

# <markdowncell>

# Verify

# <codecell>

from daltools import one

# <codecell>

Sbl = one.read(filename='S/hf_S.AOONEINT').unpack()
S = Sbl.unblock()
#print Sbl

# <markdowncell>

# Should be 16

# <codecell>

print (Da&S) + (Db&S)

# <codecell>

fa, fb = fockab(Da, Db, 'x', filename=ao2soint)

# <codecell>

#print fa - fb

# <markdowncell>

# This convention illustrates that singlet combination is zero. The fa and fb matrix should be seen as single-contraction of twoso.

# <codecell>

fabs =  [fockab(Da, Db, c, filename=ao2soint) for c in "xyz"]

# <codecell>

ls2 = [0.5*(fa - fb) for fa, fb in fabs]

# <markdowncell>

# In principle these should be added as is to one-electron

# <codecell>

oneso = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')
ls1 = prop.read(*oneso, filename=aoproper)

# <codecell>

ls = [one + two for one, two in zip(ls1, ls2)]

# <codecell>

#print ls

# <codecell>

plsp = [prefactor*p.T*M*p for M in ls]
print ls[0].shape

# <codecell>

V = makeV(plsp)

# <codecell>

#print V

# <codecell>

eigvals = numpy.linalg.eigvals(V)

# <codecell>

er = numpy.sort(eigvals.real)
print er

# <codecell>

split = abs(.5*er[:2].sum() - .25*er[2:].sum())

# <codecell>

print split*eV

# <markdowncell>

# A function for these steps

# <codecell>

def two_p_split(cmo, occa, occb):
    occupied_alpha = cmo.get_columns(occa)
    occupied_beta = cmo.get_columns(occb)
    
    Da = (occupied_alpha*occupied_alpha.T()).unblock()
    Db = (occupied_beta*occupied_alpha.T()).unblock()
    
    ls1 = prop.read(*oneso, filename=aoproper)
    ls2 = [0.5*(fa - fb) for fa, fb in [fockab(Da, Db, c, filename=ao2soint) for c in "xyz"]]
    
    ls = [one + two for one, two in zip(ls1, ls2)]
    
    plsp = [prefactor*p.T*ls_component*p for ls_component in ls]
    
    V = makeV(plsp)
    
    eigvals = numpy.linalg.eigvals(V)
    er = numpy.sort(eigvals.real)
    
    print er
    split = abs(.5*er[:2].sum() - .25*er[2:].sum())
    
    return split*eV

    

# <codecell>

two_p_split(cmo, nocc_alpha, nocc_beta)

# <codecell>

two_p_split(cmorel, nocc_alpha, nocc_beta)

# <codecell>

(1.1309485892685953-1.1287511649406907)/1.1287511649406907

# <codecell>


