#!/bin/bash -l
#
# allocate nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name 
#PBS -N 3_FBD_skyline
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME 
# change to the working directory
cd /home/woody/gwpa/gwpa007h/BEAST2_contraband/

# run your job
# java -jar contraband.jar -statefile TAXA_88_PCs_88/independent_run_3/FBD_skyline.state TAXA_88_PCs_88/independent_run_3/FBD_skyline.xml

# resume your job
java -jar contraband.jar -statefile TAXA_88_PCs_88/independent_run_3/FBD_skyline.state -resume TAXA_88_PCs_88/independent_run_3/FBD_skyline.xml
