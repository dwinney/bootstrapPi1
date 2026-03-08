#! /usr/bin/bash

if [ $# -lt 2 ]; then 
    echo "No arguments passed!"
    exit
fi

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
    # cp $BOOTSTRAP_SRC/bash/run.sh $BOOTSTRAP/fit_$i
    (
        echo '#! /usr/bin/bash '
        echo
        echo '#SBATCH --job-name=pi1_'$i
        echo '#SBATCH --output=output_'$i'.txt'
        echo '#SBATCH --partition=icn'
        echo '#SBATCH --ntasks=1'
        echo '#SBATCH --nodes=1'
        echo
        echo source /home/icn/dwinney/.bashrc
        echo module purge
        echo module load lamod/boost/1.80
        echo module load lamod/ROOT/v6-28-04
        echo srun ./do_fit
    ) > $BOOTSTRAP/fit_$i/run.sh
    chmod u+x $BOOTSTRAP/fit_$i/run.sh
    # increment
    (( i += 1 ))
done
