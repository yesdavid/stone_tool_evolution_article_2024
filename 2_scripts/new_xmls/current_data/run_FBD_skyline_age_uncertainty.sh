#!/bin/bash -l
#
# allocate nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name 
#PBS -N FBD_skyline_age_uncertainty
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME 
# change to the working directory
cd /home/woody/gwpa/gwpa007h/BEAST2_contraband/

# run your job
java -jar contraband.jar -statefile FBD_skyline_age_uncertainty.state FBD_skyline_age_uncertainty.xml

# resume your job
# java -jar contraband.jar -statefile FBD_skyline_age_uncertainty.state -resume FBD_skyline_age_uncertainty.xml
