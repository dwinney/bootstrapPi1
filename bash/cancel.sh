#! /usr/bin/bash

min=$1
max=$2

i=$min
while [[ $i -le $max ]] ; do
	scancel --name pi1_$i
	(( i += 1 ))	
done
