#! /usr/bin/bash

min=$1
max=$2

i=$min
while [[ $i -le $max ]] ; do
	rm -rf $BOOTSTRAP/fit_$i
	(( i += 1 ))
done

$BOOTSTRAP_SRC/bash/setup.sh $min $max
