#' generate_slurm_file
#'
#' @param bash_out_dir
#' @param nreps
#' @param nperms
#' @param mc.cores
#' @param J
#' @param n_j
#' @param a
#' @param b
#' @param c_p
#' @param theta_ab
#' @param optimizer
#' @param job_name
#' @param job_time
#' @param job_mem
#' @param partition
#' @param email_address
#' @param save_dir
#'
#' @return
#' @export
#'
#' @examples
generate_slurm_file <- function(bash_out_dir, nreps, nperms, mc.cores, J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2, optimizer = "bobyqa", job_name = "permediatr", job_time = "2-00:00:00", job_mem = "5G", partition = "ncf_holy", email_address = NULL, save_dir = "./"){

  path_to_script <- system.file(file.path('bin', 'permediatr_simulation.R'), package = 'permediatr')

  array_df <- expand.grid(a = a, b = b, c_p = c_p, theta_ab = theta_ab)
  job_numbers <- as.numeric(rownames(array_df))-1
  job_array_range <- paste(range(job_numbers), collapse = '-')
  ndigits <- max(floor(log10(abs(job_numbers)))+1)
  sim_names <- sprintf(paste0('pmjob_%0', ndigits,'d'), job_numbers)
  name_values <- paste(sim_names, collapse = ' ')
  param_values <- sapply(array_df, paste, collapse = ' ')

  mail_line <- NULL
  if(!is.null(email_address)){
    mail_line <- paste0('
#SBATCH --mail-user=', email_address, '
#SBATCH --mail-type=ALL')
  }

  bash_code <- paste0(
    '#!/bin/bash
#
#SBATCH -J ', job_name,'
#SBATCH --time=', job_time,'
#SBATCH -n 1
#SBATCH --cpus-per-task=',mc.cores,'
#SBATCH --mem=', job_mem,'
#SBATCH --partition=', partition,'
# Outputs ----------------------------------
#SBATCH -o %x-%A-%a.out
#SBATCH -e %x-%A-%a.err', mail_line,'
# ------------------------------------------
#
#usage: sbatch --array=', job_array_range,' ', path_to_script,'

nreps=', nreps,'
nperms=', nperms,'
mccores=', mc.cores,'
J=', J,'
nj=', n_j,'
a=(', param_values[['a']],')
b=(', param_values[['b']],')
cp=(', param_values[['c_p']],')
thetaab=', param_values[['theta_ab']],'
optimizer="', optimizer,'"
dataout="', save_dir,'"
name=(', name_values,')
index=${SLURM_ARRAY_TASK_ID}
permediatrpath=', path_to_script,'

Rscript "${permediatrpath}" --nreps "${nreps}" --nperms "${nperms}" --mc.cores "${mccores}" --J "${J}" --n_j "${nj}" --a "${a[${index}]}" --b "${b[${index}]}" --c_p "${cp[${index}]}" --theta_ab "${thetaab[${index}]}" --optimizer "${optimizer}" "${dataout}" "${name[${index}]}"')

  writeLines(bash_code, con = file.path(bash_out_dir, 'slurmjob.bash'))
}
