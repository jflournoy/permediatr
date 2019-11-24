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



