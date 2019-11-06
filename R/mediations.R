#' indirect_within.lme4
#'
#' Variables should be centered appropriately (i.e., within-person)
#'
#' @param data Data frame. Should not have missing values.
#' @param indices For permuting. Should be the same length as dim(data)[1]
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
#'
#' @return indirect and total effect.
#' @import lme4
#' @export
indirect_within.lme4 <- function(data, indices.y = NULL, indices.m = NULL, y.name, x.name, m_b.name, m_a.name, group.id, covariates.y=NULL, covariates.m=NULL, random.a=T, random.b=T, random.c_p=T, optimizer = "nloptwrap", lmeropts = list()) {
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
  starFormula.y <- update(form.y, as.formula(paste0('y_star ~ .')))
  starFormula.m <- update(form.m, as.formula(paste0('m_star ~ .')))
  maxfun <- permediatr::getmaxfun(form.y, data)
  lmeropts <- modifyList(list(na.action='na.fail',
                              REML=F,
                              control = lme4::lmerControl(optimizer = optimizer, optCtrl = list(maxfun = maxfun))),
                         lmeropts)
  if(is.null(indices.y) && is.null(indices.m)){
    #this is not a permutation
    e <- try({
      model.y <- do.call(lme4::lmer,
                         c(list(formula = form.y,
                              data = data),
                           lmeropts))

      model.m <- do.call(lme4::lmer,
                         c(list(formula = form.m,
                              data = data),
                           lmeropts))
      list(model.y, model.m)
    })
  } else if(!is.null(indices.y) && !is.null(indices.m)){
    e <- try({
      residsModel.y <- do.call(lme4::lmer,
                               c(list(formula = form.y,
                                    data = data),
                                 lmeropts))
      epsilon_z.y <- residuals(residsModel.y)
      P_j.epsilon_z.y <- epsilon_z.y[indices.y]
      Zy <- predict(residsModel.y)
      data$y_star <- P_j.epsilon_z.y + Zy
      model.y <- do.call(lme4::lmer,
                         c(list(formula = starFormula.y,
                              data = data),
                           lmeropts))
      residsModel.m <- do.call(lme4::lmer,
                               c(list(formula = form.m,
                                    data = data),
                                 lmeropts))
      epsilon_z.m <- residuals(residsModel.m)
      P_j.epsilon_z.m <- epsilon_z.m[indices.m]
      Zm <- predict(residsModel.m)
      data$m_star <- P_j.epsilon_z.m + Zm
      model.m <- do.call(lme4::lmer,
                         c(list(formula = starFormula.m,
                              data = data),
                           lmeropts))
      list(model.y, residsModel.y, model.m, residsModel.m)
    })
  } else {
    stop('Indices must be all NULL or given.')
  }
  if(inherits(e, 'try-error')){
    within.indirect.effect <- NA
    warnings <- NULL
  } else {
    if(random.a==T){
      a <- coef(model.m)[[group.id]][[x.name]]
    } else {
      a <- fixef(model.m)[x.name]
    }
    if(random.b==T){
      b <- coef(model.y)[[group.id]][[m_b.name]]
    } else {
      b <- fixef(model.y)[m_b.name]
    }
    #Within-Group Indirect Effects
    within.indirect.effect <- mean(a * b, na.rm=T)
    warnings <- unlist(lapply(e, function(fit) fit@optinfo$warnings))
  }
  return(list(ab = within.indirect.effect, warnings = warnings))
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
#'
#' @examples
permute_within <- function(n, data, group.id, series = F){
  blocks <- as.data.frame(data)[, group.id]
  requireNamespace('permute', quietly = TRUE)
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

#' getmaxfun
#'
#' @param form
#' @param data
#'
#' @return
#' @export
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
