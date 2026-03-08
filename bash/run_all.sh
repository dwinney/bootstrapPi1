#! /usr/bin/bash

# Range of folders to create
min=$1
max=$2

i=$min
while [[ $i -le $max ]] ; do
    # Go to folder
    cd $BOOTSTRAP/fit_$1
    # add to queue
    squeue batch.sh
    # increment
    (( i += 1 ))
done