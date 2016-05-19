#!/bin/bash
w=.01
s=.0001
python ~/dev/py/spectrum/gaussian.py $(cat hf_cc-pVDZ.hso1.eV) --common-width $w --step $s --output  hf_cc-pVDZ.dat1
python ~/dev/py/spectrum/gaussian.py $(cat hf_cc-pVDZ.hso2.eV) --common-width $w --step $s --output  hf_cc-pVDZ.dat2
python ~/dev/py/spectrum/gaussian.py $(cat hf_cc-pVDZ.hso.eV)  --common-width $w --step $s --output  hf_cc-pVDZ.dat

justplot2.py hf*.dat*
