#!/bin/bash
#
#SBATCH -J bootmediatr
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
#usage: sbatch --array=0-29 ../inst/bin//slurmjob.bash

#--MODULES--#
module load gcc/8.2.0-fasrc01
module load R/3.5.1-fasrc01
#-----------#

nreps=500
niter=2500
mccores=24
simtype=bootstrap
J=100
nj=10
a=(0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0)
b=(0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8 0 0 0 0.2 0.8)
cp=(0 0 0 0 0 0.8 0.8 0.8 0.8 0.8 0 0 0 0 0 0.8 0.8 0.8 0.8 0.8 0 0 0 0 0 0.8 0.8 0.8 0.8 0.8)
thetaab=(0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2)
reform=(NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA)
permtype=(groups groups groups groups groups groups groups groups groups groups within within within within within within within within within within both both both both both both both both both both)
boottype=(cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases cases)
optimizer="bobyqa"
dataout="/net/holynfs01/srv/export/mclaughlin/share_root/users/jflournoy/data/boot"
name=(boot_cases_00 boot_cases_01 boot_cases_02 boot_cases_03 boot_cases_04 boot_cases_05 boot_cases_06 boot_cases_07 boot_cases_08 boot_cases_09 boot_cases_10 boot_cases_11 boot_cases_12 boot_cases_13 boot_cases_14 boot_cases_15 boot_cases_16 boot_cases_17 boot_cases_18 boot_cases_19 boot_cases_20 boot_cases_21 boot_cases_22 boot_cases_23 boot_cases_24 boot_cases_25 boot_cases_26 boot_cases_27 boot_cases_28 boot_cases_29)
index=${SLURM_ARRAY_TASK_ID}
permediatrpath=/home/jflournoy/Rlibs/permediatr/bin/permediatr_simulation.R

Rscript "${permediatrpath}" --nreps "${nreps}" --niter "${niter}" --mc.cores "${mccores}" --simtype "${simtype}" --J "${J}" --n_j "${nj}" --a "${a[${index}]}" --b "${b[${index}]}" --c_p "${cp[${index}]}" --theta_ab "${thetaab[${index}]}" --optimizer "${optimizer}" --reform "${reform[${index}]}" --permtype "${permtype[${index}]}" --boottype "${boottype[${index}]}" "${dataout}" "${name[${index}]}"
