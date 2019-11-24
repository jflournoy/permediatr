library(argparse)
library(permediatr)

parser <- ArgumentParser(description='Run a simulation study using permediatr')
parser$add_argument('save_dir', type="character",
                    help='This is the path where the results will be saved')
parser$add_argument('simulation_name', type="character", help='Name of the simulation')
parser$add_argument('--nreps', type="integer",
                    help = 'Number of repetitions. New data will generated for each repetition.',
                    required = TRUE)
parser$add_argument('--nperms', type="integer",
                    help = 'Number of iterations. This many permutations or bootstrap simulations will be run for each repetition.')
parser$add_argument('--niter', type="integer",
                    help = 'Alternative form for "--nperms"')
parser$add_argument('--mc.cores', type="integer",
                    help = 'Number of cores to parallelize over.',
                    required = TRUE)
parser$add_argument('--simtype', type="character",
                    help = 'Either "permutation" or "bootstrap".',
                    required = TRUE)
parser$add_argument('--J', type="integer",
                    help = 'Number of groups.',
                    required = TRUE)
parser$add_argument('--n_j', type="integer",
                    help = 'Number of observations per group.',
                    required = TRUE)
parser$add_argument('--a', type="double",
                    help = 'Coefficient for `a` path (the regression of m on x).',
                    required = TRUE)
parser$add_argument('--b', type="double",
                    help = 'Coefficient for `b` path (the regression of y on m).',
                    required = TRUE)
parser$add_argument('--c_p', type="double",
                    help = 'Coefficient for `c\'` path (the regression of y on x controlling for m).',
                    required = TRUE)
parser$add_argument('--theta_ab', type="double",
                    help = 'Controls the correlation of the random effects of `a` and `b`.',
                    required = TRUE)
parser$add_argument('--optimizer', type="character", help='Name of the optimizer to use in the lme4::lmer calls.', default = 'bobyqa')
parser$add_argument('--reform', type="character", help='Formula specifying random effects structure for prediction. NULL and NA have special meaning (see lme4 documentation).', default = 'NULL')
parser$add_argument('--permtype', type="character", help='How to permute values. Valid options are "between", "within" and "between_within". "within" will shuffle values within levels of the grouping factor. "between" shuffles strata of observations defined by the grouping factor. "between_within" combines the two (note that observations are not mixed across strata).', default = 'within')
parser$add_argument('--boottype', type="character", help='What kind of bootstrapping to be done. Semiparametric is still experimental and does not work with all values of "--reform". See documentaion for lme4::bootMer.', default = 'parametric')
args <- parser$parse_args()

if(args$reform == 'NULL'){
  re.form <- NULL
} else if(args$reform == 'NA'){
  re.form <- NA
} else {
  re.form <- try(as.formula(args$reform))
}

if(inherits(re.form, 'try-error')){
  stop('Problem with --reform value: ', args$reform)
}

if(!args$permtype %in% c('between', 'within', 'between_within')){
  stop('--permtype specified incorrectly. Should be "between", "within" or "between_within", not: ', args$permtype)
}

if(is.null(args$niter) && is.null(args$nperms)){
  stop('Specify number of iterations using "--niter".')
} else if(!is.null(args$niter) && !is.null(args$nperms)) {
  stop('Specify only one of "--nperms" or "--niter"')
} else {
  if(is.null(args$niter)){
    niter <- args$niter <- args$nperms
  } else {
    niter <- args$nperms <- args$niter
  }
}

results <- permediatr::run_simulation(nreps = as.numeric(args$nreps),
                                      niter = niter,
                                      mc.cores = args$mc.cores,
                                      simtype = args$simtype,
                                      J = args$J,
                                      n_j = args$n_j,
                                      a = args$a,
                                      b = args$b,
                                      c_p = args$c_p,
                                      theta_ab = args$theta_ab,
                                      optimizer = args$optimizer,
                                      re.form = re.form,
                                      permtype = args$permtype,
                                      boottype = args$boottype)

rds_file <- file.path(args$save_dir, paste0(args$simulation_name, '.rds'))
message('Saving results to RDS file: ', rds_file)
saveRDS(results, file = rds_file)
csv_file <- file.path(args$save_dir, paste0(args$simulation_name, '.csv'))
message('Saving results to CSV file: ', csv_file)
results_df <- do.call(rbind,lapply(results, function(x) data.frame(t(unlist(x[-which(names(x) %in% c('ab_warnings', 'ab_singular'))])))))
results_csv <- cbind(results_df, data.frame(args))
write.csv(results_csv, file = csv_file)

# args <- parser$parse_args(
#   c('~/', 'test',
#     '--nreps', '2',
#   '--niter', '500',
#   '--mc.cores', '7',
#   '--simtype', 'bootstrap',
#   '--J', '100',
#   '--n_j', '10',
#   '--a', '.2',
#   '--b', '0',
#   '--c_p', '0',
#   '--theta_ab', '.2',
#   '--optimizer', 'bobyqa',
#   '--reform', 'NA',
#   '--permtype', 'between_within',
#   '--boottype', 'parametric')
# )
