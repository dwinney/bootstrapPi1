#! /usr/bin/bash

# Range of folders to create
min=$1
max=$2

i=$min
while [[ $i -le $max ]] ; do
    # Create new fit directory
    mkdir $BOOTSTRAP/fit_$i
    # Copy a version of the executable
    cp $BOOTSTRAP_SRC/bin/do_fit   $BOOTSTRAP/fit_$i
    # Also copy a version of the batch script
    cp $BOOTSTRAP_SRC/bash/batch.sh $BOOTSTRAP/fit_$i
    # increment
    (( i += 1 ))
done