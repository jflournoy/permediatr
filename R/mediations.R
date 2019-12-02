#' lmemm_within
#'
#' @param data
#' @param y.name y variable name.
#' @param x.name x variable name (should be within-person centered).
#' @param m_b.name mediator variable name (within-person centered for computing b path).
#' @param m_a.name mediator variable name (raw for computing a path).
#' @param group.id grouping variable for random effects.
#' @param covariates.y any covariates to include in regression of y.
#' @param covariates.m any covariates to include in regression of m.
#' @param random.a include random effect for m on x (a path).
#' @param random.b include random effect for y on m (b path).
#' @param random.c_p include random effect for y on x (c' path).
#' @param optimizer
#' @param lmeropts
#'
#' @return
#' @export
#' @import lme4
#' @examples
lmemm_within <- function(data, y.name, x.name, m_b.name, m_a.name, group.id, covariates.y=NULL, covariates.m=NULL, random.a=T, random.b=T, random.c_p=T, optimizer = "nloptwrap", lmeropts = list()){
  requireNamespace('lme4', quietly = TRUE)
  re_terms.y <- c('1', c(m_b.name, x.name)[c(random.b, random.c_p)])
  re_terms.m <- c('1', c(x.name)[random.a])
  re_part.y <- paste0('(', paste(re_terms.y, collapse = ' + '), ' | ', group.id, ')')
  re_part.m <- paste0('(', paste(re_terms.m, collapse = ' + '), ' | ', group.id, ')')
  form.y <- as.formula(paste0(y.name, ' ~ 1 + ',
                              paste(covariates.y, collapse = ' + '),
                              ' + ', m_b.name, ' + ', x.name,
                              ' + ', re_part.y))
  form.m <- as.formula(paste0(m_a.name, ' ~ 1 + ',
                              paste(covariates.m, collapse = ' + '),
                              ' + ', x.name,
                              ' + ', re_part.m))
  maxfun <- permediatr::getmaxfun(form.y, data)
  lmeropts <- modifyList(list(na.action='na.fail',
                              REML=F,
                              control = lme4::lmerControl(optimizer = optimizer, optCtrl = list(maxfun = maxfun))),
                         lmeropts)
  e <- try({
    model.y <- do.call(lme4::lmer,
                       c(list(formula = form.y,
                              data = data),
                         lmeropts))

    model.m <- do.call(lme4::lmer,
                       c(list(formula = form.m,
                              data = data),
                         lmeropts))
    r <- list(model.y = model.y,
              model.m = model.m,
              lmeropts = lmeropts,
              data = data,
              parspec = c(y.name = y.name,
                          x.name = x.name,
                          m_b.name = m_b.name,
                          m_a.name = m_a.name,
                          group.id = group.id,
                          covariates.y = covariates.y,
                          covariates.m = covariates.m,
                          random.a = random.a,
                          random.b = random.b,
                          random.c_p = random.c_p,
                          optimizer = optimizer))
    attr(r,'class') <- 'permediatrMod'
    r
  })
  return(e)
}

#' permediatrPermMod
#'
#' @param x a permediatrMod object
#' @param indices.y indices to use to permute residuals for model.y
#' @param indices.m indices to use to permute residuals for model.m
#' @param re.form.y random effects form for y
#' @param re.form.m random effects form for m
#'
#' @return an object of class permediatrPermMod
#' @export
#' @import lme4
#' @examples
permediatrPermMod <- function(x, indices.y = NULL, indices.m = NULL, re.form.y = NULL, re.form.m = NULL){
  form.y <- formula(x$model.y)
  form.m <- formula(x$model.m)

  starFormula.y <- update(form.y, as.formula(paste0('y_star ~ .')))
  starFormula.m <- update(form.m, as.formula(paste0('m_star ~ .')))

  data <- x$data

  e <- try({
    residsModel.y <- x$model.y
    Zy <- predict(residsModel.y, re.form = re.form.y)
    epsilon_z.y <- residsModel.y@frame$y - Zy
    P_j.epsilon_z.y <- epsilon_z.y[indices.y]
    data$y_star <- P_j.epsilon_z.y + Zy
    model.y <- do.call(lme4::lmer,
                       c(list(formula = starFormula.y,
                              data = data),
                         x$lmeropts))
    residsModel.m <- x$model.m
    Zm <- predict(residsModel.m, re.form = re.form.m)
    epsilon_z.m <- residsModel.m@frame$m - Zm
    P_j.epsilon_z.m <- epsilon_z.m[indices.m]
    data$m_star <- P_j.epsilon_z.m + Zm
    model.m <- do.call(lme4::lmer,
                       c(list(formula = starFormula.m,
                              data = data),
                         x$lmeropts))
    r <- list(model.y = model.y,
              model.m = model.m,
              residsModel.y = residsModel.y,
              residsModel.m = residsModel.m,
              lmeropts = x$lmeropts,
              data = data,
              parspec = c(x$parspec,
                          indices.y = indices.y,
                          indices.m = indices.m,
                          re.form.y = re.form.y,
                          re.form.m = re.form.m))
    attr(r,'class') <- c('permediatrMod', 'permediatrPermMod')
    r
  })
  return(e)
}

#' is.permediatrMod
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.permediatrMod <- function(x) {
  inherits(x, 'permediatrMod')
}

#' is.permediatrPermMod
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.permediatrPermMod <- function(x) {
  inherits(x, 'permediatrPermMod')
}

#' indirectEffect
#'
#' Takes a permediatrMod or permediatrPermMod object and extracts the estimate
#' of the indirect effect.
#'
#' @param x a permediatrMod object
#'
#' @return an indirectEffect object, with elements ab, the indirect effect
#'   estiamte, warnings, a vector of warnings, and singular, a vector of
#'   logicals indicating if convergence was singular.
#' @export
#' @import lme4
#'
#' @examples
indirectEffect <- function(x){
  requireNamespace('lme4', quietly = TRUE)
  if(inherits(x, 'try-error')){
    within.indirect.effect <- NA
    warnings <- NULL
    singulars <- NULL
    ab_obj <- list(ab = within.indirect.effect, warnings = warnings, singular = singulars)
  } else if(is.permediatrMod(x)) {
    ab_obj <- with(as.list(x$parspec), {
      if(random.a==T){
        a <- coef(x$model.m)[[group.id]][[x.name]]
      } else {
        a <- fixef(x$model.m)[x.name]
      }
      if(random.b==T){
        b <- coef(x$model.y)[[group.id]][[m_b.name]]
      } else {
        b <- fixef(x$model.y)[m_b.name]
      }
      #Within-Group Indirect Effects
      within.indirect.effect <- mean(a * b, na.rm=T)

      if(is.permediatrPermMod(x)){
        model_list <- x[c('model.y', 'model.m', 'residsModel.y', 'residsModel.m')]
      } else {
        model_list <- x[c('model.y', 'model.m')]
      }
      singulars <- unlist(lapply(model_list, function(fit) {
        s <- NULL
        if(isSingular(fit)){
          s <- TRUE
        }
        return(s)
      }))
      warnings <- unlist(lapply(model_list, function(fit) {
        warnings <- fit@optinfo$warnings
        return(warnings)
      }))
      list(ab = within.indirect.effect, warnings = warnings, singular = singulars)
    })
  } else {
    stop('Variable "x" has incorrect class: ', class(x))
  }
  attr(ab_obj,'class') <- 'indirectEffect'
  return(ab_obj)
}

#' lmemm_within.perm
#'
#' Variables should be centered appropriately (i.e., within-person)
#'
#'
#' @param ... Passed on to the permediatrPermMod constructor
#'
#' @return indirectEffect object.
#' @export
lmemm_within.perm <- function(...) {
  ab_perm_mod <- permediatrPermMod(...)
  ab <- indirectEffect(ab_perm_mod)
  return(ab)
}

#' permute_within
#'
#' @param n
#' @param data
#' @param group.id
#' @param series
#'
#' @return
#' @export
#' @import permute
#'
#' @examples
permute_within <- function(n, data, group.id, series = F){
  requireNamespace('permute', quietly = TRUE)
  blocks <- as.data.frame(data)[, group.id]
  if(series){
    type <- 'series'
  } else{
    type <- 'free'
  }
  n_obs <- dim(data)[1]
  perm_control <- permute::how(nperm = n, within = permute::Within(type = type), blocks = blocks)
  p <- permute::shuffleSet(n = n_obs, control = perm_control)
  return(p)
}

#' permute_between
#'
#' @param n
#' @param data
#' @param group.id
#' @param series
#'
#' @return
#' @export
#' @import permute
#'
#' @examples
permute_between <- function(n, data, group.id, series = F){
  requireNamespace('permute', quietly = TRUE)
  strata <- as.data.frame(data)[, group.id]
  n_obs <- dim(data)[1]
  perm_control <- permute::how(nperm = n, within = permute::Within(type = 'none'), plots = permute::Plots(strata = strata, type = 'free'), blocks = NULL)
  p <- permute::shuffleSet(n = n_obs, control = perm_control)
  return(p)
}

#' permute_between_within
#'
#' @param n
#' @param data
#' @param group.id
#' @param series
#'
#' @return
#' @export
#' @import permute
#'
#' @examples
permute_between_within <- function(n, data, group.id, series = F){
  requireNamespace('permute', quietly = TRUE)
  strata <- as.data.frame(data)[, group.id]
  if(series){
    type <- 'series'
  } else{
    type <- 'free'
  }
  n_obs <- dim(data)[1]
  perm_control <- permute::how(nperm = n, within = permute::Within(type = type), plots = permute::Plots(strata = strata, type = 'free'), blocks = NULL)
  p <- permute::shuffleSet(n = n_obs, control = perm_control)
  return(p)
}

#' getmaxfun
#'
#' @param form
#' @param data
#'
#' @return
#' @export
#' @import lme4
#'
#' @examples
getmaxfun <- function(form, data){
  requireNamespace('lme4', quietly = TRUE)
  model_form <- lme4::lFormula(form, data = data)
  model_numFx <- length(dimnames(model_form$X)[[2]])
  model_numRx <- sum(as.numeric(lapply(model_form$reTrms$cnms, function(x) {
    l <- length(x)
    (l*(l - 1)) / 2 + l
  })))
  model_maxfun <- 10*(model_numFx + model_numRx + 1)^2
  return(model_maxfun)
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

  ab_bc <- permediatr:::bias_corrected_est(est = ab$ab, ab_vec)

  dists <- list(c_ab = ab_vec - ab$ab,
                c_mean = ab_vec - ab_bc,
                neg_c_ab = -(ab_vec - ab$ab),
                neg_c_mean = -(ab_vec - ab_bc))
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

#' bias_corrected_est
#'
#' @param est
#' @param replicates
#' @param na.rm
#'
#' @return numeric bias corrected estimate.
#'
bias_corrected_est <- function(est, replicates, na.rm = T){
  bc <- 2*est - mean(replicates, na.rm = na.rm)
  return(bc)
}

p_from_replicates <- function(est, replicates){
  est_bc <- permediatr:::bias_corrected_est(est = est, replicates = replicates)
  dists <- list(c_est = replicates - est,
                c_estbc = replicates - est_bc,
                neg_c_est = -(replicates - est),
                neg_c_estbc = -(replicates - est_bc))

  p_dists_est <- lapply(1:length(dists), function(i){
    p <- try(ecdf(dists[[i]])(est))
    if(inherits(p, 'try-error')){
      warning(p)
      warning('Error in computing p value: ', names(dists)[[i]], '\nest: ', est)
      p <- NA
    }
    return(p)
  })
  names(p_dists_est) <- paste0('p_est_', names(dists))
  p_dists_estbc <- lapply(1:length(dists), function(i){
    p <- try(ecdf(dists[[i]])(est_bc))
    if(inherits(p, 'try-error')){
      warning(p)
      warning('Error in computing p value: ', names(dists)[[i]], '\nest: ', est)
      p <- NA
    }
    return(p)
  })
  names(p_dists_estbc) <- paste0('p_estbc_', names(dists))

  p_ci <- try(min(sum(replicates < 0), sum(replicates > 0))/length(replicates))
  if(inherits(p_ci, 'try-error')){
    warning(p_ci)
    warning('Error in computing p value: p_ci')
    p_ci <- NA
  }
  r <- c(p_dists_est, p_dists_estbc, list(p_ci = p_ci))
  return(r)
}

#' bootMerCases
#'
#' @param x a fitted \code{merMod} object: see \code{\link[lme4]{lmer}},
#'   \code{\link[lme4]{glmer}}, etc.
#' @param FUN a function taking a fitted \code{merMod} object as input and
#'   returning the statistic of interest, which must be a (possibly named)
#'   numeric vector.
#' @param nsim number of resampled data sets, positive integer; the bootstrap B
#'   (or R).
#' @param seed optional argument to \code{\link[base]{set.seed}}.
#' @param how Character specifying whether resampling should be done of groups,
#'   observations within groups, or both. Currently only one level of grouping
#'   is allowed.
#' @param parallel The type of parallel operation to be used (if any).
#' @param ncpus integer: number of processes to be used in parallel operation:
#'   typically one would choose this to be the number of available CPUs.
#'
#' @return
#' @import data.table
#' @export
#'
bootMerCases <- function(x, FUN, nsim = 1, seed = NULL, how = c('groups', 'within', 'both'), parallel = c("no", "multicore"), ncpus = getOption("boot.ncpus", 1L)){
  requireNamespace('data.table', quietly = TRUE)
  data.table::setDTthreads(ncpus)
  how <- grep(how, c('groups', 'within', 'both'), value = T)
  parallel <- grep(parallel, c("no", "multicore"), value = T)

  if(!how %in% c('groups', 'within', 'both'))
    stop("'how' not set correctly.")
  model_data <- x@frame
  clusters <- c(rev(names(lme4::getME(x, "flist"))), ".id")

  a_resample <- function(data, clusters, how){
    data <- as.data.table(data)
    index <- 1:dim(data)[1]
    group_ids <- data[, get(clusters[1])]
    unique_ids <- unique(group_ids)
    if(how %in% c('groups', 'both')){
      resampled_ids <- data.table::setnames(data.table(sample(unique_ids, size = length(unique_ids), replace = TRUE)),
                                            clusters[1])
      resampled_ids[, new_id := 1:.N]
      resampled_data <- data[resampled_ids, on = clusters[1]]
      data.table::setnames(resampled_data,
                           old = c('new_id', clusters[1]),
                           new = c(clusters[1], 'old_id'))[, old_id := NULL]
    } else {
      resampled_data <- data
    }
    if(how %in% c('within', 'both')){
      newrows <- resampled_data[, sample(.I, size = .N, replace = TRUE), by = eval(clusters[1])][[2]]
      resampled_data <- resampled_data[newrows]
    }
    return(as.data.frame(resampled_data))
  }

  #the data sets are generated serially, so we don't have to worry about
  #multicore leading to irreproducible results.
  set.seed(seed)
  rep.data <- lapply(integer(nsim), function(sim) a_resample(data = model_data, clusters = clusters, how = how))

  t0 <- FUN(x)

  if(parallel == 'multicore'){
    split_rep.data <- suppressWarnings(split(rep.data, 1:mc.cores))
    tstar <- unlist(parallel::mclapply(split_rep.data, function(D_list) {
      lapply(D_list, function(D){
        FUN(update(x, data = D))
      })
    }, mc.cores = ncpus), recursive = FALSE)
  } else {
    tstar <- lapply(rep.data, function(D) {
      FUN(update(x, data = D))
    })
  }

  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)

  RES <- structure(list(t0 = t0, t = t(tstar), R = nsim, data = model_data,
                        seed = .Random.seed, statistic = FUN,
                        sim = "case", call = match.call()),
                   class = "boot")

  return(RES)
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
#'   Not used if type == 'cases'.
#' @param type character string specifying the type of bootstrap,
#'   \code{"parametric"} or \code{"semiparametric"}; partial matching is
#'   allowed.
#' @param how character specifying the type of resampling if type == 'cases'.
#' @param seed
#'
#' @return
#' @export
#' @import lme4
#'
#' @examples
boot_ab <- function(x, nsim, ncpus, re.form = NA, type = c('parametric', 'semiparametric', 'cases'), how = c('groups', 'within', 'both'), seed = NULL){

  # Where random-number generation is done in the worker processes, the default
  # behaviour is that each worker chooses a separate seed, non-reproducibly.
  # However, with parallel = "multicore" or parallel = "snow" using the default
  # cluster, a second approach is used if RNGkind("L'Ecuyer-CMRG") has been
  # selected. In that approach each worker gets a different subsequence of the
  # RNG stream based on the seed at the time the worker is spawned and so the
  # results will be reproducible if ncpus is unchanged, and for parallel =
  # "multicore" if parallel::mc.reset.stream() is called: see the examples for
  # mclapply.
  requireNamespace('lme4', quietly = TRUE)
  ab <- indirectEffect(x)

  afun.a <- permediatr::make_bootFUN(x, path = 'a')
  afun.b <- permediatr::make_bootFUN(x, path = 'b')

  if(type %in% c('parametric', 'semiparametric')){
    aboot.a <- lme4::bootMer(x$model.m, FUN = afun.a, nsim = nsim, re.form = re.form, type = type, seed = seed, parallel = 'multicore', ncpus = ncpus)
    aboot.b <- lme4::bootMer(x$model.y, FUN = afun.b, nsim = nsim, re.form = re.form, type = type, seed = seed, parallel = 'multicore', ncpus = ncpus)
  } else if (type == 'cases'){
    aboot.a <- permediatr::bootMerCases(x$model.m, FUN = afun.a, nsim = nsim, seed = seed, how = how, parallel = 'multicore', ncpus = ncpus)
    aboot.b <- permediatr::bootMerCases(x$model.y, FUN = afun.b, nsim = nsim, seed = seed, how = how, parallel = 'multicore', ncpus = ncpus)
  }

  t_list.a <- lapply(seq_len(nrow(aboot.a$t)), function(i) aboot.a$t[i,])
  t_list.b <- lapply(seq_len(nrow(aboot.b$t)), function(i) aboot.b$t[i,])

  ab_vec <- mapply(function(a, b){
    mean(a*b, na.rm=T)
  }, t_list.a, t_list.b)

  r <- list(model.m = aboot.a, model.y = aboot.b, ab_t = ab_vec, ab = ab)
  attr(r, 'class') <- 'permediatrBoot'
  return(r)
}

#' is.permediatrBoot
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.permediatrBoot <- function(x) {
  inherits(x, 'permediatrBoot')
}

#' summary.permediatrBoot
#'
#' @param x An object of class permediatrBoot.
#' @param which
#' @param verbose
#'
#' @return some p-values and other information.
#' @export
#'
summary.permediatrBoot <- function(x, which = 'ab', verbose = F){
  model.m <- x$model.m
  model.y <- x$model.y
  ab_t <- x$ab_t
  ab <- x$ab
  msgs <- list(msg.a = attr(model.m, 'boot.all.msgs'),
               msg.b = attr(model.y, 'boot.all.msgs'))
  if(verbose){
    msg.rpt <- lapply(msgs, function(x) lapply(x, function(y) if(!length(y) == 0) message(names(y))))
  }
  n_conv_warn <- sum(unlist(lapply(msgs, function(x) sum(as.numeric(x[['factory-warning']])))))
  n_singular <- sum(unlist(lapply(msgs, function(x) sum(as.numeric(x[['factory-message']])))))
  message('Number of convergence warnings: ', n_conv_warn, '\n')
  message('Number of singular convergence: ', n_singular, '\n')

  if(which == 'ab'){
    est <- ab$ab
    replicates <- ab_t
  } else if(which == 'a'){
    est <- mean(model.m$t0)
    replicates <- apply(model.m$t, 1, mean)
  } else if(which == 'b'){
    est <- mean(model.y$t0)
    replicates <- apply(model.y$t, 1, mean)
  }
  est_bc <- permediatr:::bias_corrected_est(est = est, replicates = replicates)
  p_values <- p_from_replicates(est = est, replicates = replicates)
  results_list <- c(
    parname = which,
    p_values,
    list(
      n_conv_warn = n_conv_warn,
      n_singular = n_singular,
      est = est,
      est_bc = est_bc,
      warnings = ab$warnings,
      singular = ab$singular,
      mean_t = mean(replicates)),
    as.list(quantile(replicates, probs = c(.025, .975, seq(0, 1, 0.25)))))
  return(results_list)
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


