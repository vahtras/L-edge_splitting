#!/bin/bash
w=${width-.01}
s=${step-.0001}

for i in $@
do
    gaussian.py $(cat $i) --common-width $w --step $s --output $i.dat
done

