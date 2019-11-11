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
                    help = 'Number of permutations. This many permutations will be run for each repetition.',
                    required = TRUE)
parser$add_argument('--mc.cores', type="integer",
                    help = 'Number of cores to parallelize over.',
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
args <- parser$parse_args()

results <- run_permutation_simulation(nreps = as.numeric(args$nreps),
                                      nperms = args$nperms,
                                      mc.cores = args$mc.cores,
                                      J = args$J,
                                      n_j = args$n_j,
                                      a = args$a,
                                      b = args$b,
                                      c_p = args$c_p,
                                      theta_ab = args$theta_ab,
                                      optimizer = args$optimizer)
rds_file <- file.path(args$save_dir, paste0(args$simulation_name, '.rds'))
message('Saving results to RDS file: ', rds_file)
saveRDS(results, file = rds_file)
csv_file <- file.path(args$save_dir, paste0(args$simulation_name, '.csv'))
message('Saving results to CSV file: ', csv_file)
results_df <- do.call(rbind,lapply(results, function(x) data.frame(t(unlist(x[-which(names(x) == 'ab_warnings')])))))
results_csv <- cbind(results_df, data.frame(args))
write.csv(results_csv, file = csv_file)
