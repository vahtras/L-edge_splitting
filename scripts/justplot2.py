#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as mpl
import sys


STYLES = {'solid': '-', 'dash': '--', 'dot': ':'}


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='Files with data')
    parser.add_argument('--title', help='Title')
    parser.add_argument('--suptitle', help='Suptitle')
    parser.add_argument('--xlabel', help='xlabel')
    parser.add_argument('--ylabel', help='ylabel')
    parser.add_argument('--filenames', action='store_true', help='filenames')
    parser.add_argument('--legends', nargs='*', help='Filenames or custom legends')
    parser.add_argument('--colors', help='one-character colors [bgrcmyk]')
    parser.add_argument('--styles', nargs='+',choices=['solid', 'dash', 'dot'], help='Line styles')
    parser.add_argument('--xlim', nargs=2, type=float, help='x-axis range')
    parser.add_argument('--loc', help='Legend location')


    args = parser.parse_args()

    if args.legends:
        assert len(args.files) == len(args.legends), "Different number of files/legends"
        legends = args.legends
    elif args.filenames:
        legends = args.filenames
    else:
        legends = ' '*len(args.files)

    if args.colors:
        assert len(args.files) == len(args.colors), "Different number of files/colors"
        colors = args.colors
    else:
        colors = 'bgrcmyk'

    if args.styles:
        assert len(args.files) == len(args.styles), "Different number of files/styles"
        styles = args.styles
    else:
        styles = ('solid',)*len(args.files)

    if args.suptitle:
        mpl.suptitle(args.suptitle)

    if args.title:
        mpl.title(args.title)

    if args.xlabel:
        mpl.xlabel(args.xlabel, fontsize=14)

    if args.ylabel:
        mpl.ylabel(args.ylabel, fontsize=14)

    if args.xlim:
        mpl.xlim(*args.xlim)


    for f,l,col,s in zip(args.files, legends, colors, styles):
        npf = np.loadtxt(f)
        for c in range(1, npf.shape[1]):
            mpl.plot(npf[:, 0], npf[:, c], label=l, color=col, linestyle=STYLES[s])

    if args.filenames or args.legends:
        if args.loc:
            mpl.legend(loc=args.loc)
        else:
            mpl.legend()

    mpl.show()

if __name__ == "__main__":
    sys.exit(main())
