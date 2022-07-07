#!/bin/bash -l
#
# allocate 4 nodes
#PBS -l nodes=1:ppn=4,walltime=WALLTIME_PLACEHOLDER
#
# job name 
#PBS -N JOBNAME_PLACEHOLDER
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME 
# change to the working directory
cd /home/woody/gwpa/gwpa007h/BEAST2_contraband/

# run your job
java -jar contraband.jar SCRIPTPATH_PLACEHOLDER
