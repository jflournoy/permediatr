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
#' @param nreps
#' @param mc.cores
#' @param J
#' @param n_j
#' @param a
#' @param b
#' @param c_p
#' @param theta_ab
#' @param optimizer
#'
#' @return
#' @export
#' @import pbapply
#'
#' @examples
run_permutation_simulation <- function(nreps, nperms, mc.cores, J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2, optimizer = "bobyqa"){
  message("\nRunning ", nreps, " simulations with the following parameters:
permutations: ", nperms,"
J: ", J,"
n_j: ", n_j,"
a: ", a,"
b: ", b,"
c_p: ", c_p,"
theta_ab: ", theta_ab,"
optimizer: ", optimizer)
  #could make this slightly more efficient by splitting up reps differently.
  nreps_digits <- 1 + floor(log10(nreps))
  pbapply::pboptions(type = 'timer', char = '+', style = 3)
  reps <- pbapply::pblapply(1:nreps, function (i) {
    message(sprintf(paste0('\nReplication % ', nreps_digits, 'd out of %d'), i, nreps))
    message('Generating new data...')
    adf <- simulate_mediation_data(J = J, n_j = n_j, a = a, b = b, c_p = c_p, theta_ab = theta_ab)
    ab <- indirect_within.lme4(data = adf,
                               indices.y = NULL,
                               indices.m = NULL,
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
    message('Generating permutations...')
    perms.y <- permute_within(n = nperms, data = adf, group.id = 'id', series = F)
    perms.m <- permute_within(n = nperms, data = adf, group.id = 'id', series = F)
    perms_list.y <- lapply(seq_len(nrow(perms.y)), function(i) perms.y[i,])
    perms_list.m <- lapply(seq_len(nrow(perms.m)), function(i) perms.m[i,])
    split_perms_list.y <- suppressWarnings(split(perms_list.y, 1:mc.cores))
    split_perms_list.m <- suppressWarnings(split(perms_list.m, 1:mc.cores))
    message('Evaluating ', nperms, ' permutations over ', mc.cores, ' processors...')
    perm_set_time <- system.time({
      ab_perms <- parallel::mcmapply(function(perm_list.y, perm_list.m){
        mapply(function(indices.y, indices.m){
          ab <- indirect_within.lme4(data = adf,
                                     indices.y = indices.y,
                                     indices.m = indices.m,
                                     y.name = 'y',
                                     x.name = 'x_wc',
                                     m_b.name = 'm_wc',
                                     m_a.name = 'm',
                                     group.id = 'id',
                                     covariates.y=c('cov.y', 'm_bc', 'x_bc'),
                                     covariates.m=c('cov.m', 'x_bc'),
                                     random.a=T,
                                     random.b=T,
                                     random.c_p=T,
                                     optimizer = optimizer)
        }, perm_list.y, perm_list.m, SIMPLIFY = FALSE)
      }, split_perms_list.y, split_perms_list.m, SIMPLIFY = FALSE, mc.cores = mc.cores)
    })
    message('Parallel evaluitation took: ', perm_set_time[['elapsed']], 's')
    no_convwarnings <- unlist(lapply(unlist(ab_perms, recursive = F), function(aperm){
      awarning <- aperm[['warnings']]
      is.null(awarning)
    }))
    ab_vec <- unlist(lapply(unlist(ab_perms, recursive = F), `[[`, 'ab'))
    ab_vec <- ab_vec[no_convwarnings]
    ab_vec_centered_mean <- ab_vec - mean(ab_vec)
    ab_vec_centered_ab <- ab_vec - ab$ab
    p_ab_perm_mean_c <- ecdf(ab_vec_centered_mean)(ab$ab)
    p_ab_perm_ab_c <- ecdf(ab_vec_centered_ab)(ab$ab)
    p_opsign_ci <- mean(sign(ab$ab)*ab_vec < 0)
    n_conv_warn <- sum(!no_convwarnings)
    message('Finished this replication. Number of convergence warnings: ', n_conv_warn, '\n')
    results_list <- list(p_ab_perm_mean_c = p_ab_perm_mean_c,
         p_ab_perm_ab_c = p_ab_perm_ab_c,
         p_opsign_ci = p_opsign_ci,
         n_conv_warn = n_conv_warn,
         ab = ab$ab,
         ab_warnings = ab$warnings,
         mean_ab_perm = mean(ab_vec))
    return(results_list)
  })
  return(reps)
}
