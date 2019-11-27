#' simulate_mediation_data
#'
#' @param J
#' @param n_j
#' @param a
#' @param b
#' @param c_p
#' @param theta_ab
#'
#' @return
#' @export
#' @import lme4
#'
#' @examples
simulate_mediation_data <- function(J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2){
  requireNamespace('lme4', quietly = TRUE)
  re_a = T
  re_b = T
  newdata <- data.frame(x = rnorm(J*n_j),
                        cov.m = rep(rnorm(J), each = n_j),
                        cov.y = rep(rnorm(J), each = n_j),
                        id = rep(factor(1:J), each = n_j))
  if(re_a){
    re_term.m <- '(1 + x | id)'
    if(re_b){
      re_term.y <- '(1 + m + x | id)'
    } else {
      re_term.y <- re_term.m
    }
  } else {
    re_term.m <- '(1 | id)'
    if(re_b){
      re_term.y <- '(1 + m | id)'
    } else {
      re_term.y <- re_term.m
    }
  }
  form.m <- as.formula(paste0('~ 1 + cov.m + x + ', re_term.m))
  params.m <- list(theta = c(`id.(Intercept)` = 1,
                             `id.x.(Intercept)` = .1,
                             `id.x` = .5),
                   beta = c(`(Intercept)` = 0,
                            cov.m = 0,
                            x = a),
                   sigma = 1)
  newdata$m <- simulate(form.m,
                        newdata = newdata,
                        newparams = params.m,
                        family = 'gaussian', re.form = NA)[[1]]
  form.y <- as.formula(paste0('~ 1 + cov.y + m + x + ', re_term.y))
  params.y <- list(theta = c(`id.(Intercept)` = 1,
                             `id.m.(Intercept)` = .1,
                             `id.x.(Intercept)` = .1,
                             `id.m` = .5,
                             `id.x.m` = theta_ab,
                             `id.x` = .5),
                   beta = c(`(Intercept)` = 0,
                            cov.y = 0,
                            x = c_p,
                            m = b),
                   sigma = 1)
  newdata$y <- simulate(form.y,
                        newdata = newdata,
                        newparams = params.y,
                        family = 'gaussian', re.form = NA)[[1]]
  #Center data within person, and center within-person means
  x_w_mean <- aggregate(newdata$x, by = list(newdata$id), FUN = mean, na.rm = T)
  colnames(x_w_mean) <- c("id", "x_w_mean")
  newdata <- merge(newdata, x_w_mean, by="id")
  newdata$x_mean <- mean(x_w_mean$x_w_mean)
  newdata$x_bc <- newdata$x_w_mean - newdata$x_mean
  newdata$x_wc <- newdata$x - newdata$x_w_mean
  m_w_mean <- aggregate(newdata$m, by = list(newdata$id), FUN = mean, na.rm = T)
  colnames(m_w_mean) <- c("id", "m_w_mean")
  newdata <- merge(newdata, m_w_mean, by="id")
  newdata$m_mean <- mean(m_w_mean$m_w_mean)
  newdata$m_bc <- newdata$m_w_mean - newdata$m_mean
  newdata$m_wc <- newdata$m - newdata$m_w_mean
  return(newdata)
}

#' run_permutation_simulation
#'
#' Included for backward compatibility.
#'
#' @param nreps Number of repetitions using new data.
#' @param niter For each simulated data set, number of permutation or bootstrap iterations.
#' @param mc.cores Number of cores to use.
#' @param ... To be passed to \code{\link{run_simulation}}
#'
#' @return
#' @export
#'
#' @examples
run_permutation_simulation <- function(nreps, niter, mc.cores, ...){
  return(run_simulation(nreps = nreps, niter = niter, mc.cores = mc.cores, simtype = 'permutation'))
}



#' run_simulation
#'
#' @param nreps Number of repetitions using new data.
#' @param niter For each simulated data set, number of permutation or bootstrap
#'   iterations.
#' @param mc.cores Number of cores to use.
#' @param simtype Permutations or bootstrap.
#' @param J
#' @param n_j
#' @param a
#' @param b
#' @param c_p
#' @param theta_ab
#' @param optimizer
#' @param re.form
#' @param permtype This can be used to specify either the permutation type
#'   ('within', 'between', or 'between_within'), or if simtype is 'bootstrap',
#'   and boottype is 'cases', it is used to specify how cases are resampled,
#'   ('groups', 'within', 'both').
#' @param boottype
#'
#' @return
#' @export
#' @import pbapply
#'
#' @examples
run_simulation <- function(nreps, niter, mc.cores, simtype = 'permutation', J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2, optimizer = "bobyqa", re.form = NULL, permtype = 'within', boottype = 'parametric'){
  simtype <- grep(simtype, c('permutation', 'bootstrap'), value = TRUE)
  thecall <- match.call()
  message("\nRunning ", nreps, " simulations with the following parameters:
", paste(mapply(paste, names(thecall), thecall, MoreArgs = list(sep = ': ')), collapse = '\n'))

  #could make this slightly more efficient by splitting up reps differently.
  nreps_digits <- 1 + floor(log10(nreps))
  pbapply::pboptions(type = 'timer', char = '+', style = 3)
  reps <- pbapply::pblapply(1:nreps, function (i) {
    message(sprintf(paste0('\nReplication % ', nreps_digits, 'd out of %d'), i, nreps))
    message('Generating new data...')
    adf <- simulate_mediation_data(J = J, n_j = n_j, a = a, b = b, c_p = c_p, theta_ab = theta_ab)
    ab_mod <- permediatr::lmemm_within(data = adf,
                                       y.name = 'y',
                                       x.name = 'x_wc',
                                       m_b.name = 'm_wc',
                                       m_a.name = 'm',
                                       group.id = 'id',
                                       covariates.y=c('cov.y', 'm_bc', 'x_bc'),
                                       covariates.m=c('cov.m', 'x_bc'),
                                       random.a=T,
                                       random.b=T,
                                       random.c_p=T, optimizer = optimizer)
    if(simtype == 'permutation'){
      results_list <- permute_ab(x = ab_mod, nperms = niter, mc.cores = mc.cores, re.form = re.form)
    } else if(simtype == 'bootstrap'){
      theboot <- boot_ab(x = ab_mod, nsim = niter, ncpus = mc.cores, re.form = re.form, type = boottype, how = permtype)
      results_list <- do.call(rbind,
                              list(summary(theboot, which = 'ab', verbose = T),
                                   summary(theboot, which = 'a'),
                                   summary(theboot, which = 'b')))
    }
    return(results_list)
  })
  return(reps)
}


#' generate_slurm_file
#'
#' @param bash_out_dir
#' @param nreps
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
#' @param niter
#' @param simtype
#' @param re.form
#'
#' @return
#' @export
#'
#' @examples
generate_slurm_file <- function(bash_out_dir, nreps, niter, mc.cores, simtype = 'permutation', J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2, re.form = NULL, permtype = 'within', boottype = 'parametric', optimizer = "bobyqa", job_name = "permediatr", job_time = "2-00:00:00", job_mem = "5G", partition = "ncf_holy", email_address = NULL, save_dir = "./", fname_prefix = 'pmjob'){

  path_to_script <- system.file(file.path('bin', 'permediatr_simulation.R'), package = 'permediatr')

  array_df <- expand.grid(a = a, b = b, c_p = c_p, theta_ab = theta_ab, re.form = re.form, permtype = permtype, boottype = boottype)
  array_df <- array_df[array_df$a == 0 | array_df$b == 0, ]
  job_numbers <- 1:dim(array_df)[1]-1
  job_array_range <- paste(range(job_numbers), collapse = '-')
  ndigits <- max(floor(log10(abs(job_numbers)))+1)
  sim_names <- sprintf(paste0(fname_prefix, '_%0', ndigits,'d'), job_numbers)
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
#usage: sbatch --array=', job_array_range,' ', file.path(bash_out_dir, 'slurmjob.bash'),'

#--MODULES--#
module load gcc/8.2.0-fasrc01
module load R/3.5.1-fasrc01
#-----------#

nreps=', nreps,'
niter=', niter,'
mccores=', mc.cores,'
simtype=', simtype,'
J=', J,'
nj=', n_j,'
a=(', param_values[['a']],')
b=(', param_values[['b']],')
cp=(', param_values[['c_p']],')
thetaab=(', param_values[['theta_ab']],')
reform=(', param_values[['re.form']],')
permtype=(', param_values[['permtype']],')
boottype=(', param_values[['boottype']],')
optimizer="', optimizer,'"
dataout="', save_dir,'"
name=(', name_values,')
index=${SLURM_ARRAY_TASK_ID}
permediatrpath=', path_to_script,'

Rscript "${permediatrpath}" --nreps "${nreps}" --niter "${niter}" --mc.cores "${mccores}" --simtype "${simtype}" --J "${J}" --n_j "${nj}" --a "${a[${index}]}" --b "${b[${index}]}" --c_p "${cp[${index}]}" --theta_ab "${thetaab[${index}]}" --optimizer "${optimizer}" --reform "${reform[${index}]}" --permtype "${permtype[${index}]}" --boottype "${boottype[${index}]}" "${dataout}" "${name[${index}]}"')

  writeLines(bash_code, con = file.path(bash_out_dir, 'slurmjob.bash'))
}
