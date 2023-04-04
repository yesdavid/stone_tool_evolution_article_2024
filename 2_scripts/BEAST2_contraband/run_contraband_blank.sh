#!/bin/bash -l
#
# allocate nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=WALLTIME_PLACEHOLDER
#
# job name 
#PBS -N JOBNAME_PLACEHOLDER
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME 
# change to the working directory
cd /home/woody/gwpa/gwpa007h/BEAST2_contraband/

# run your job
java -jar contraband.jar -statefile SCRIPTPATH_PLACEHOLDER.state SCRIPTPATH_PLACEHOLDER

# resume your job
# java -jar contraband.jar -statefile SCRIPTPATH_PLACEHOLDER.state -resume SCRIPTPATH_PLACEHOLDER
