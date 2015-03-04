from dalmisc.twoso import *
from util import full
D = full.matrix((34, 34))
D[0, 0] = 2
F = fock(D, component='x', filename='hf_S.AO2SOINT')

print F
print F[31, 5]

Da = full.matrix((34, 34))
Db = Da
Da[0, 0] =  1
fa, fb = fockab(Da, Db, component='x', filename='hf_S.AO2SOINT')
print .5*(fa-fb)


