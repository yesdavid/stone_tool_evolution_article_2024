#!/bin/bash -l
#
# allocate 4 nodes
#PBS -l nodes=1:ppn=4,walltime=24:00:00
#
# job name 
#PBS -N BMPruneLikelihood_calval_1_FBDbds_BDMMprime_TAXA71_PCs12_nIter1000m
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME 
# change to the working directory
cd /home/woody/gwpa/gwpa007h/BEAST2_contraband/

# run your job
java -jar contraband.jar ./TAXA71_PCs12/script_BMPruneLikelihood_calval_1_FBDbds_BDMMprime_TAXA71_PCs12_nIter1000m.xml
