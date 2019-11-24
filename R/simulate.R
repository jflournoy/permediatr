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

#' permute_ab
#'
#' @param x a permediatrMod object
#' @param nperms
#' @param mc.cores
#' @param re.form
#' @param permtype
#'
#' @return
#' @export
#'
#' @examples
permute_ab <- function(x, nperms, mc.cores, re.form = NULL, permtype = 'within'){
  ab <- permediatr::indirectEffect(x)
  adf <- x$data
  group.id <- x$parspec[['group.id']]
  permute_fun <- get(paste0('permute_', permtype), envir = as.environment('package:permediatr'))
  message('Generating permutations...')
  perms.y <- permute_fun(n = nperms, data = adf, group.id = group.id, series = F)
  perms.m <- permute_fun(n = nperms, data = adf, group.id = group.id, series = F)
  perms_list.y <- lapply(seq_len(nrow(perms.y)), function(i) perms.y[i,])
  perms_list.m <- lapply(seq_len(nrow(perms.m)), function(i) perms.m[i,])
  split_perms_list.y <- suppressWarnings(split(perms_list.y, 1:mc.cores))
  split_perms_list.m <- suppressWarnings(split(perms_list.m, 1:mc.cores))
  message('Evaluating ', nperms, ' permutations over ', mc.cores, ' processors...')
  perm_set_time <- system.time({
    ab_perms <- parallel::mcmapply(function(perm_list.y, perm_list.m){
      mapply(function(indices.y, indices.m){
        ab <- permediatr::lmemm_within.perm(x = x,
                                            indices.y = indices.y,
                                            indices.m = indices.m,
                                            re.form.y = re.form,
                                            re.form.m = re.form)
      }, perm_list.y, perm_list.m, SIMPLIFY = FALSE)
    }, split_perms_list.y, split_perms_list.m, SIMPLIFY = FALSE, mc.cores = mc.cores)
  })
  message('Parallel evaluitation took: ', perm_set_time[['elapsed']], 's')
  no_convwarnings <- unlist(lapply(unlist(ab_perms, recursive = F), function(aperm){
    awarning <- aperm[['warnings']]
    is.null(awarning)
  }))
  no_singular <- unlist(lapply(unlist(ab_perms, recursive = F), function(aperm){
    s <- aperm[['singular']]
    is.null(s)
  }))
  ab_vec <- unlist(lapply(unlist(ab_perms, recursive = F), `[[`, 'ab'))
  ab_vec <- ab_vec[no_convwarnings]

  dists <- list(c_ab = ab_vec - ab$ab,
                c_mean = ab_vec - mean(ab_vec),
                neg_c_ab = -(ab_vec - ab$ab),
                neg_c_mean = -(ab_vec - mean(ab_vec)))
  p_dists <- lapply(1:length(dists), function(i){
    p <- try(ecdf(dists[[i]])(ab$ab))
    if(inherits(p, 'try-error')){
      warning(p)
      warning('Error in computing p value: ', names(dists)[[i]])
      p <- NA
    }
    return(p)
  })
  names(p_dists) <- paste0('p_', names(dists))

  p_opsign_ci <- try(mean(sign(ab$ab)*ab_vec < 0))
  if(inherits(p_opsign_ci, 'try-error')){
    warning(p_opsign_ci)
    warning('Error in computing p value: p_opsign_ci')
    p_opsign_ci <- NA
  }

  n_conv_warn <- sum(!no_convwarnings)
  n_singular <- sum(!no_singular)
  message('Finished this set of permutations. Number of convergence warnings: ', n_conv_warn, '\n')
  message('Finished this set of permutations. Number of singular convergence: ', n_singular, '\n')
  results_list <- c(
    p_dists,
    list(
      p_opsign_ci = p_opsign_ci,
      n_conv_warn = n_conv_warn,
      n_singular = n_singular,
      ab = ab$ab,
      ab_warnings = ab$warnings,
      ab_singular = ab$singular,
      mean_ab_perm = mean(ab_vec)))
  return(results_list)
}

#' boot_ab
#'
#' See \code{\link[lme4]{bootMer}} for more information on arguments for
#' \code{re.form} and \code{type}.
#'
#' @param x a permediatrMod object.
#' @param nsim
#' @param ncpus
#' @param re.form formula, NA (equivalent to \code{\link[lme4]{bootMer}}'s
#'   \code{use.u=FALSE}; new normal deviates are drawn), or NULL (equivalent to
#'   \code{\link[lme4]{bootMer}}'s \code{use.u=TRUE}; all inference is
#'   conditional on these values) or formula specifying which random effects to
#'   incorporate. See \code{\link[lme4]{simulate.merMod}} for further details.
#' @param type character string specifying the type of bootstrap,
#'   \code{"parametric"} or \code{"semiparametric"}; partial matching is
#'   allowed.
#'
#' @return
#' @export
#' @import lme4
#'
#' @examples
boot_ab <- function(x, nsim, ncpus, re.form = NA, type = c('parametric', 'semiparametric')){
  requireNamespace('lme4', quietly = TRUE)
  ab <- indirectEffect(x)

  afun.a <- permediatr::make_bootFUN(x, path = 'a')
  afun.b <- permediatr::make_bootFUN(x, path = 'b')

  aboot.a <- lme4::bootMer(x$model.m, FUN = afun.a, nsim = nsim, re.form = re.form, type = type, parallel = 'multicore', ncpus = ncpus)
  aboot.b <- lme4::bootMer(x$model.y, FUN = afun.b, nsim = nsim, re.form = re.form, type = type, parallel = 'multicore', ncpus = ncpus)

  t_list.a <- lapply(seq_len(nrow(aboot.a$t)), function(i) aboot.a$t[i,])
  t_list.b <- lapply(seq_len(nrow(aboot.b$t)), function(i) aboot.b$t[i,])

  ab_vec <- mapply(function(a, b){
    mean(a*b, na.rm=T)
  }, t_list.a, t_list.b)

  msgs <- list(msg.a = attr(aboot.a, 'boot.all.msgs'),
               msg.b = attr(aboot.b, 'boot.all.msgs'))
  msg.rpt <- lapply(msgs, function(x) lapply(x, function(y) if(!length(y) == 0) message(names(y))))

  n_conv_warn <- sum(unlist(lapply(msgs, function(x) sum(as.numeric(x[['factory-warning']])))))
  n_singular <- sum(unlist(lapply(msgs, function(x) sum(as.numeric(x[['factory-message']])))))
  message('Finished this set of bootstrapping. Number of convergence warnings: ', n_conv_warn, '\n')
  message('Finished this set of bootstrapping. Number of singular convergence: ', n_singular, '\n')

  boot_ab_est <- mean(ab_vec)

  dists <- list(c_ab = ab_vec - ab$ab,
                c_mean = ab_vec - boot_ab_est,
                neg_c_ab = -(ab_vec - ab$ab),
                neg_c_mean = -(ab_vec - boot_ab_est))
  p_dists_ab <- lapply(1:length(dists), function(i){
    p <- try(ecdf(dists[[i]])(ab$ab))
    if(inherits(p, 'try-error')){
      warning(p)
      warning('Error in computing ab p value: ', names(dists)[[i]])
      p <- NA
    }
    return(p)
  })
  names(p_dists_ab) <- paste0('p_ab_', names(dists))
  p_dists_boot <- lapply(1:length(dists), function(i){
    p <- try(ecdf(dists[[i]])(boot_ab_est))
    if(inherits(p, 'try-error')){
      warning(p)
      warning('Error in computing boot\'d ab p value: ', names(dists)[[i]])
      p <- NA
    }
    return(p)
  })
  names(p_dists_boot) <- paste0('p_boot_', names(dists))

  p_opsign_ci <- try(mean(sign(boot_ab_est)*ab_vec < 0))
  if(inherits(p_opsign_ci, 'try-error')){
    warning(p_opsign_ci)
    warning('Error in computing p value: p_opsign_ci')
    p_opsign_ci <- NA
  }
  results_list <- c(
    p_dists_ab,
    p_dists_boot,
    list(
      p_opsign_ci = p_opsign_ci,
      n_conv_warn = n_conv_warn,
      n_singular = n_singular,
      ab = ab$ab,
      boot_ab = boot_ab_est,
      ab_warnings = ab$warnings,
      ab_singular = ab$singular),
    as.list(quantile(ab_vec, probs = c(.025, .975, seq(0, 1, 0.25)))))
  return(results_list)
}

#' run_simulation
#'
#' @param nreps Number of repetitions using new data.
#' @param niter For each simulated data set, number of permutation or bootstrap iterations.
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
#' @param permtype
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
      results_list <- boot_ab(x = ab_mod, nsim = niter, ncpus = mc.cores, re.form = re.form, type = boottype)
    }
    return(results_list)
  })
  return(reps)
}

#' make_bootFUN
#'
#' @param x A permediatrMod object with the full mediation specification.
#' @param path Which path is being bootstrapped, the 'a' path, or the 'b' path?
#' @param is.null.model Should the boot FUN expect to get a null model without
#'   the effect of interest? Set this to true of you're doing parametric or
#'   semiparametric bootstrapping for the nil null distribution. The
#'   \code{\link{make_null_formula}} should be used to create the null model to
#'   be the basis for (semi-)parametric simulation.
#'
#' @return
#' @export
#'
#' @examples
make_bootFUN <- function(x, path = c('a', 'b'), is.null.model = F){
  func <- with(as.list(x$parspec), {
    if(path == 'a'){
      name <- x.name
      random <- random.a
      o.name <- m_a.name
      #In the case that the model passed to the bootstrap FUN is a null model
      #based on bootstrap simulations, we want to use the formula from the
      #original model and estimate the effect of interest in the context of that
      #null model. The og model will be passed to the updator and updated if
      #is.null.model is true.
      og <- x$model.m
    } else if(path == 'b'){
      name <- m_b.name
      random <- random.b
      o.name <- y.name
      og <- x$model.y
    } else {
      stop('Incorrect path specification')
    }
    func <- function(model){
      if(is.null.model){
        newdata__ <- og@frame
        newdata__[, o.name] <- model@frame[,o.name]
        umodel <- update(model, formula(og), data = newdata__)
      } else {
        umodel <- model
      }
      if(random==T){
        r <- coef(umodel)[[group.id]][[name]]

      } else {
        r <- fixef(umodel)[name]
      }
      return(r)
    }
    func
  })
  return(func)
}

#' make_null_model
#'
#' @param x A permediatrMod object with the full mediation specification.
#' @param path For which path is the null model being generated, the 'a' path, or the 'b' path?
#'
#' @return
#' @export
#'
#' @examples
make_null_model <- function(x, path = c('a', 'b')){
  mod <- with(as.list(x$parspec), {
    if(path == 'a'){
      name <- x.name
      random <- random.a
      o.name <- m_a.name
      og <- x$model.m
    } else if(path == 'b'){
      name <- m_b.name
      random <- random.b
      o.name <- y.name
      og <- x$model.y
    } else {
      stop('Incorrect path specification')
    }
    mod <- update(og, permediatr::make_null_formula(model = og,
                                                    name = name,
                                                    random = random))
    mod
  })
  return(mod)
}

#' make_null_formula
#'
#' @param model A specific model from a permediatrMod object.
#' @param name The variable name to be excluded from the model.
#' @param random Whether the variable name is also included in the random
#'   effects part of the model.
#'
#' @return
#' @export
#'
#' @examples
make_null_formula <- function(model, name, random){
  fixed_form <- nobars(formula(model))
  if(random){
    re_list <- findbars(formula(model))
    re_null <- paste(lapply(re_list, function(x) {
      re <- as.character(re_list[[1]])
      re_2 <- gsub(paste0('(.*?) \\+ ', name,'(.*?)'), '\\1\\2', re[[2]])
      paste0('(', re_2, ' ', re[[1]], ' ', re[[3]],')')
    }), collapse = ' + ')
    form_null <- update(fixed_form, paste(c(paste0('. ~ . -', name), re_null), collapse = ' + '))
  } else {
    form_null <- update(formula(model), paste0('. ~ . -', name))
  }
  return(form_null)
}
