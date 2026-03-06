#! /usr/bin/bash

# Range of folders to create
min=1
max=10

i=$min
while [[ $i -le $max ]] ; do
    mkdir fit_$i
    (( i += 1 ))
done