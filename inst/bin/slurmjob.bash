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
#usage: sbatch --array=0-9 ../inst/bin//slurmjob.bash

#--MODULES--#
module load gcc/8.2.0-fasrc01
module load R/3.5.1-fasrc01
#-----------#

nreps=500
nperms=2500
mccores=24
simtype=bootstrap
J=100
nj=10
a=(0 0.2 0.8 0 0 0 0.2 0.8 0 0)
b=(0 0 0 0.2 0.8 0 0 0 0.2 0.8)
cp=(0 0 0 0 0 0.8 0.8 0.8 0.8 0.8)
thetaab=(0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2)
reform=(NA NA NA NA NA NA NA NA NA NA)
permtype=(within within within within within within within within within within)
boottype=(parametric parametric parametric parametric parametric parametric parametric parametric parametric parametric)
optimizer="bobyqa"
dataout="/net/holynfs01/srv/export/mclaughlin/share_root/users/jflournoy/data"
name=(pmjob_0 pmjob_1 pmjob_2 pmjob_3 pmjob_4 pmjob_5 pmjob_6 pmjob_7 pmjob_8 pmjob_9)
index=${SLURM_ARRAY_TASK_ID}
permediatrpath=/home/jflournoy/Rlibs/permediatr/bin/permediatr_simulation.R

Rscript "${permediatrpath}" --nreps "${nreps}" --niter "${niter}" --mc.cores "${mccores}" --simtype "${simtype}" --J "${J}" --n_j "${nj}" --a "${a[${index}]}" --b "${b[${index}]}" --c_p "${cp[${index}]}" --theta_ab "${thetaab[${index}]}" --optimizer "${optimizer}" --reform "${reform[${index}]}" --permtype "${permtype[${index}]}" --boottype "${boottype[${index}]}" "${dataout}" "${name[${index}]}"
