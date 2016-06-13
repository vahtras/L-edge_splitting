#!/usr/bin/env python
import sys
import numpy
from .constants import CM, EV

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
        conv = EV
    elif args.cm:
        unit  = "cm"
        conv = CM 
    else:
        unit = "a.u"
        conv = 1

    print((args.fmt + " %s") % (mean(args.n1, args.n2, args.filename)*conv, unit))
         

if __name__ == "__main__":
    main()
    
