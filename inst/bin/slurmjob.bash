#!/bin/bash
#
#SBATCH -J permediatr
#SBATCH --time=3-00:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH --mem=10G
#SBATCH --partition=ncf_holy
# Outputs ----------------------------------
#SBATCH -o %x-%A-%a.out
#SBATCH -e %x-%A-%a.err
#SBATCH --mail-user=john_flournoy@fas.harvard.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
#
#usage: sbatch --array=0-33 ../inst/bin//slurmjob.bash

nreps=500
nperms=2500
mccores=24
J=100
nj=10
a=(0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0)
b=(0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8)
cp=(0 0 0 0 0 0.8 0.8 0.8 0.8 0.8 0 0 0 0 0 0.8 0.8 0.8 0.8 0.8)
thetaab=(0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2)
reform=(NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL NA NA NA NA NA NA NA NA NA NA)
optimizer="bobyqa"
dataout="/net/holynfs01/srv/export/mclaughlin/share_root/users/jflournoy/data"
name=(pmjob_00 pmjob_01 pmjob_02 pmjob_03 pmjob_06 pmjob_09 pmjob_10 pmjob_11 pmjob_12 pmjob_15 pmjob_18 pmjob_19 pmjob_20 pmjob_21 pmjob_24 pmjob_27 pmjob_28 pmjob_29 pmjob_30 pmjob_33)
index=${SLURM_ARRAY_TASK_ID}
permediatrpath=/home/jflournoy/Rlibs/permediatr/bin/permediatr_simulation.R

Rscript "${permediatrpath}" --nreps "${nreps}" --nperms "${nperms}" --mc.cores "${mccores}" --J "${J}" --n_j "${nj}" --a "${a[${index}]}" --b "${b[${index}]}" --c_p "${cp[${index}]}" --theta_ab "${thetaab[${index}]}" --optimizer "${optimizer}" --reform "${reform[${index}]}" "${dataout}" "${name[${index}]}"
