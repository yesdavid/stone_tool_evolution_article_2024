#!/bin/bash -l
#
# allocate nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name 
#PBS -N RUN_PLACEHOLDER:RUN_SCRIPT_PLACEHOLDER
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME 
# change to the working directory
# cd /home/

# run your job
# java -jar contraband.jar -statefile CURRENT_FOLDER_PLACEHOLDER/CURRENT_SCRIPT_PLACEHOLDER.state CURRENT_FOLDER_PLACEHOLDER/CURRENT_SCRIPT_PLACEHOLDER.xml

# resume your job
java -jar contraband.jar -statefile CURRENT_FOLDER_PLACEHOLDER/CURRENT_SCRIPT_PLACEHOLDER.state -resume CURRENT_FOLDER_PLACEHOLDER/CURRENT_SCRIPT_PLACEHOLDER.xml
