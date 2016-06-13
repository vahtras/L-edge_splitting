#!/usr/bin/env python
import sys
import numpy
from scipy.constants import physical_constants as pc
eV = pc['Hartree energy in eV'][0]
cm = pc['hartree-inverse meter relationship'][0]/100

def mean(n1, n2, filename):
    eigenvalues = numpy.loadtxt(filename)
    return meansplit(n1, n2, eigenvalues)

def meansplit(n1, n2, eigenvalues):
    assert n1 + n2 == eigenvalues.shape[0]
    split = abs(sum(eigenvalues[:n1])/n1 - sum(eigenvalues[n1:])/(n2))
    return split

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('n1', type=int)
    parser.add_argument('n2', type=int)
    parser.add_argument('filename')
    parser.add_argument('--eV', action='store_true')
    parser.add_argument('--cm', action='store_true')
    parser.add_argument('--fmt', default='%f')


    args = parser.parse_args()

    if args.eV:
        unit  = "eV"
        conv = eV
    elif args.cm:
        unit  = "cm"
        conv = cm 
    else:
        unit = "a.u"
        conv = 1

    print((args.fmt + " %s") % (mean(args.n1, args.n2, args.filename)*conv, unit))
         

if __name__ == "__main__":
    main()
    
