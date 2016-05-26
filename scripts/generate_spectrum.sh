#!/bin/bash
w=.01
s=.0001

for i in $@
do
    python ~/dev/py/spectrum/gaussian.py $(cat $i) --common-width $w --step $s
    --output  $i.dat
done

justplot2.py $(for i in $@; do echo $i.dat; done)
